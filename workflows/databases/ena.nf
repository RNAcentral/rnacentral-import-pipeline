process list_subdirs {
  tag { "$name" }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '1h'

  input:
  tuple val(name), val(path)

  output:
  tuple val(name), path("${name}-dirs.txt")

  when: params.databases.ena?.run

  script:
  """
  find -L ${path} -mindepth 1 -maxdepth 1 -type d > ${name}-dirs.txt
  """
}

// The cheap prefilter: a source's signature covers the names, sizes and mtimes of
// its archives, which is everything a stat gives us for free. ENA's .ncr.gz files
// are years old, so most sources come back identical and are never fetched at all.
// Batched exactly like the fetch so the NFS walk fans out instead of running serially.
process stat_sources {
  tag { "$name" }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '4h'

  input:
  tuple val(name), val(roots)

  output:
  path('sources.tsv')

  when: params.databases.ena?.run

  script:
  """
  cat > roots.txt <<'ROOTS'
${roots.join('\n')}
ROOTS

  : > sources.tsv
  while IFS= read -r root; do
    signature=\$(find -L "\$root" -type f \( -name '*.ncr.gz' -o -name '*.tar' \) \
      -printf '%P\t%s\t%T@\n' \
      | LC_ALL=C sort \
      | sha256sum \
      | cut -d' ' -f1)
    printf '%s\t%s\n' "\$root" "\$signature" >> sources.tsv
  done < roots.txt
  """
}

// Single reduce over every source signature: what to fetch, and what the record-level
// diff is allowed to retire.
process ena_source_diff {
  input:
  path('sources*.tsv')

  output:
  path('to_fetch.txt'), emit: to_fetch
  path('scanned.txt'), emit: scanned
  path('sources.csv'), emit: sources

  when: params.databases.ena?.run

  script:
  def force = params.force_full_import ? '--force-full' : ''
  """
  cat sources*.tsv > all-sources.tsv
  rnac ena source-diff $force all-sources.tsv to_fetch.txt scanned.txt sources.csv
  """
}

process fetch_directory {
  tag { "${remotes.size()} sources" }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '2d'

  input:
  val(remotes)

  output:
  path('chunks/*.ncr'), optional: true

  when: params.databases.ena?.run

  script:
  """
  cat > roots.txt <<'ROOTS'
${remotes.join('\n')}
ROOTS

  mkdir chunks

  # One source at a time, each split into chunks named after it: that name is all a
  # chunk keeps of where its records came from, and it is what lets a later run skip
  # an unchanged source without its records looking deleted. The copy is thrown away
  # between sources so one batch cannot fill the work directory. The outer loop reads
  # on fd 3, leaving stdin free for tar and zcat.
  while IFS= read -r root <&3; do
    label=\$(printf '%s' "\$root" | sha1sum | cut -c1-16)

    rm -rf copied
    rsync \
      -aL --partial \
      --prune-empty-dirs \
      --include='*/' \
      --include='**/*.ncr.gz' \
      --include='**/*.tar' \
      --exclude='*.fasta.gz' \
      "\$root" copied

    find copied -type f -empty -delete

    # A truncated archive in the snapshot must not cost the whole batch, so unpack each
    # one on its own and carry on, rewinding the .ncr over whatever a failed member
    # managed to write. Piping the lot into one xargs aborts everything with exit 123.
    skipped=0

    pushd copied
    while IFS= read -r archive; do
      if ! tar -xvf "\$archive"; then
        echo "WARN: unreadable tar, skipping \$archive" >&2
        skipped=\$(( skipped + 1 ))
      fi
    done < <(find . -type f -name '*.tar')
    popd

    : > \$label.ncr
    while IFS= read -r archive; do
      kept=\$(wc -c < \$label.ncr)
      if ! zcat "\$archive" >> \$label.ncr; then
        echo "WARN: unreadable archive, skipping \$archive" >&2
        truncate -s "\$kept" \$label.ncr
        skipped=\$(( skipped + 1 ))
      fi
    done < <(find copied -type f -name '*.ncr.gz')

    echo "\$skipped unreadable archives skipped in \$root" >&2

    if [ -s \$label.ncr ]; then
      split-ena --max-sequences ${params.databases.ena.max_sequences} \$label.ncr chunks
    else
      echo "No .ncr data fetched for \$root; emitting no chunks" >&2
    fi

    rm -rf copied \$label.ncr
  done 3< roots.txt
  """
}

process fetch_metadata {
  input:
  path(urls)

  output:
  tuple path('tpa.tsv'), path('model-lengths.csv')

  when: params.databases.ena?.run

  script:
  """
  cat $urls | xargs -I {} wget -O - {} >> tpa.tsv
  cmstat \$RIBODIR/models/ribo.0p20.extra.cm | grep -v '^#' | awk '{ printf("%s,%d\\n", \$2, \$6); }' > model-lengths.csv
  """
}

process process_file {
  memory '8GB'
  tag { "$raw" }
  time '10m'

  input:
  tuple path(raw), path(to_parse), path(tpa), path(model_lengths)

  output:
  path('*.{csv,parquet}'), optional: true

  script:
  """
  # Delta: keep only the new/changed records (KEEP_ALL on the first delta run), so
  # ribotyper and the parse run on the changed subset alone. An empty result means
  # the whole chunk was unchanged -- emit nothing. See docs/incremental-parsing-ena.md.
  kept=\$(rnac ena filter --only $to_parse $raw filtered.ncr)
  if [[ "\$kept" -eq 0 ]]; then
    echo "No new/changed records in $raw; skipping ribotyper and parse" >&2
    exit 0
  fi

  ena2fasta.py filtered.ncr sequences.fasta
  if [[ -e sequences.fasta ]]; then
    /rna/ribovore/ribotyper.pl sequences.fasta ribotyper-results
  else
    mkdir ribotyper-results
  fi
  rnac ena parse --counts $raw-counts.txt filtered.ncr $tpa ribotyper-results $model_lengths .

  mkdir $baseDir/ena-counts 2>/dev/null || true
  cp $raw-counts.txt $baseDir/ena-counts/
  """
}

// Cheap per-chunk pass: accession,signature for every record, no ribotyper or DB.
process ena_signatures {
  tag { "$raw" }
  time '30m'

  input:
  path(raw)

  output:
  path('signatures.csv')

  when: params.databases.ena?.run

  script:
  """
  rnac ena signatures $raw signatures.csv
  """
}

// Single reduce: diff all chunk signatures against the stored ENA manifest in the
// database, producing the to-parse filter plus the manifest.csv / deletions.csv the
// generic delta wiring already consumes.
process ena_delta_diff {
  memory '4GB'

  input:
  path('signatures*.csv')
  path(scanned)

  output:
  path('to_parse.txt'), emit: to_parse
  path('deletions.csv'), emit: deletions
  path('manifest.csv'), emit: manifest

  when: params.databases.ena?.run

  script:
  def force = params.force_full_import ? '--force-full' : ''
  """
  # No chunks at all when every source was unchanged; the diff still has to run, so
  # a source that has gone from the snapshot gets its records retired.
  cat signatures*.csv > all-signatures.csv 2>/dev/null || : > all-signatures.csv
  rnac ena delta-diff $force all-signatures.csv $scanned to_parse.txt deletions.csv manifest.csv
  """
}

workflow ena {
  main:
    channel.fromPath('files/import-data/ena/tpa-urls.txt') | set { urls }
    fetch_metadata(urls) | set { metadata }

    channel.fromList([
      ['wgs', "$params.databases.ena.remote/wgs/"],
      ['tls', "$params.databases.ena.remote/tls/"],
      ['tsa', "$params.databases.ena.remote/tsa/"],
    ]) \
    | list_subdirs \
    | flatMap { name, listing ->
        listing.readLines()
          .findAll { line -> line.trim() }
          .collate( params.databases.ena.subdir_batch_size )
          .collect { batch -> [name, batch.collect { s -> s.trim() }] }
      } \
    | set { subdir_batches }

    // con/ and std/ are single directories rather than trees of project subdirs, so
    // each is one source in its own right.
    channel.fromList([
      ['con', ["$params.databases.ena.remote/con/"]],
      ['std', ["$params.databases.ena.remote/std/"]],
    ]) \
    | mix( subdir_batches ) \
    | stat_sources \
    | collect \
    | set { source_signatures }

    ena_source_diff(source_signatures)

    // Only the sources whose stat signature moved are worth fetching; the rest are
    // never copied, decompressed, split or signatured. Re-batched here because the
    // changed set is normally far smaller than the listing it came from.
    ena_source_diff.out.to_fetch \
    | splitText \
    | map { line -> line.trim() } \
    | filter { line -> line } \
    | collate( params.databases.ena.subdir_batch_size ) \
    | fetch_directory \
    | flatten \
    | set { chunks }

    // Signature every chunk, reduce to one global diff. A run where nothing changed
    // produces no chunks and so no signatures, but the diff still runs.
    chunks | ena_signatures | collect | ifEmpty([]) | set { signatures }
    ena_delta_diff(signatures, ena_source_diff.out.scanned)

    // Parse only the changed records: pair every chunk with the single to-parse
    // filter and the metadata, then filter+ribotyper+parse inside process_file.
    chunks \
    | combine(ena_delta_diff.out.to_parse) \
    | combine(metadata) \
    | process_file \
    | set { parsed }

    // manifest.csv, sources.csv and deletions.csv join the data stream; import-data.nf
    // routes them (deletions -> load_deletions, the other two -> rnac manifest apply
    // once the release has committed).
    parsed \
    | mix(
      ena_delta_diff.out.manifest,
      ena_delta_diff.out.deletions,
      ena_source_diff.out.sources,
    ) \
    | set { data }

  emit: data
}
