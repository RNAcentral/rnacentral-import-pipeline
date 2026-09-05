process list_subdirs {
  tag { "$name" }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '1h'

  input:
  tuple val(name), val(path)

  output:
  tuple val(name), path("${name}-dirs.txt")

  when: params.databases.ena.run

  script:
  """
  find -L ${path} -mindepth 1 -maxdepth 1 -type d > ${name}-dirs.txt
  """
}

process fetch_directory {
  tag { "$name" }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '2d'

  input:
  tuple val(name), val(remotes)

  output:
  path("${name}-chunks/*.ncr"), optional: true

  when: params.databases.ena?.run

  script:
  """
  rsync \
    -aL --partial \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --include='**/*.tar' \
    --exclude='*.fasta.gz' \
    ${remotes.join(' ')} "copied"

  find copied -type f -empty -delete

  # A truncated archive in the snapshot must not cost the whole batch, so unpack each
  # one on its own and carry on, rewinding ${name}.ncr over whatever a failed member
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

  : > ${name}.ncr
  while IFS= read -r archive; do
    kept=\$(wc -c < ${name}.ncr)
    if ! zcat "\$archive" >> ${name}.ncr; then
      echo "WARN: unreadable archive, skipping \$archive" >&2
      truncate -s "\$kept" ${name}.ncr
      skipped=\$(( skipped + 1 ))
    fi
  done < <(find copied -type f -name '*.ncr.gz')

  echo "\$skipped unreadable archives skipped in ${name}" >&2

  mkdir $name-chunks
  if [ -s ${name}.ncr ]; then
    split-ena --max-sequences ${params.databases.ena.max_sequences} ${name}.ncr ${name}-chunks
  else
    echo "No .ncr data fetched for ${name} batch; emitting no chunks" >&2
  fi
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

  when: { params.databases.ena?.run }

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

  output:
  path('to_parse.txt'), emit: to_parse
  path('deletions.csv'), emit: deletions
  path('manifest.csv'), emit: manifest

  when: { params.databases.ena?.run }

  script:
  """
  cat signatures*.csv > all-signatures.csv
  rnac ena delta-diff all-signatures.csv to_parse.txt deletions.csv manifest.csv
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

    channel.fromList([
      ['con', ["$params.databases.ena.remote/con/"]],
      ['std', ["$params.databases.ena.remote/std/"]],
    ]) \
    | mix( subdir_batches ) \
    | fetch_directory \
    | flatten \
    | set { chunks }

    // Signature every chunk, reduce to one global diff.
    chunks | ena_signatures | collect | ena_delta_diff

    // Parse only the changed records: pair every chunk with the single to-parse
    // filter and the metadata, then filter+ribotyper+parse inside process_file.
    chunks \
    | combine(ena_delta_diff.out.to_parse) \
    | combine(metadata) \
    | process_file \
    | set { parsed }

    // manifest.csv and deletions.csv join the data stream; import-data.nf routes
    // them (deletions -> load_deletions, manifest -> rnac manifest apply post-release).
    parsed \
    | mix(ena_delta_diff.out.manifest, ena_delta_diff.out.deletions) \
    | set { data }

  emit: data
}
