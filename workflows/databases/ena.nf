process list_subdirs {
  tag { "$name" }
  when { params.databases.ena.run }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '1h'

  input:
  tuple val(name), val(path)

  output:
  tuple val(name), path("${name}-dirs.txt")

  """
  find -L ${path} -mindepth 1 -maxdepth 1 -type d -printf '%s\\t%p\\n' > ${name}-dirs.txt
  """
}

process fetch_directory {
  tag { "$name" }
  when { params.databases.ena.run }
  queue 'datamover'
  containerOptions "${params.common_container} --bind /nfs:/nfs"
  time '2d'

  input:
  tuple val(name), val(remotes)

  output:
  path("${name}-chunks/*.ncr"), optional: true

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
  find copied -type f -name '*.gz' | xargs -r -I {} gzip --quiet -l {} | awk '{ if (\$2 == 0) print \$4 }' | xargs -r -I {} rm {}.gz

  pushd copied
  find . -type f -name '*.tar' | xargs -r -I {} tar -xvf {}
  popd
  find copied -type f -name '*.ncr.gz' | xargs -r zcat > ${name}.ncr

  mkdir $name-chunks
  if [ -s ${name}.ncr ]; then
    split-ena --max-sequences ${params.databases.ena.max_sequences} ${name}.ncr ${name}-chunks
  else
    echo "No .ncr data fetched for ${name} batch; emitting no chunks" >&2
  fi
  """
}

process fetch_metadata {
  when { params.databases.ena.run }

  input:
  path(urls)

  output:
  tuple path('tpa.tsv'), path('model-lengths.csv')

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
  tuple path(raw), path(tpa), path(model_lengths)

  output:
  path('*.csv')

  """
  ena2fasta.py $raw sequences.fasta
  if [[ -e sequences.fasta ]]; then
    /rna/ribovore/ribotyper.pl sequences.fasta ribotyper-results
  else
    mkdir ribotyper-results
  fi
  rnac ena parse --counts $raw-counts.txt $raw $tpa ribotyper-results $model_lengths .

  mkdir $baseDir/ena-counts 2>/dev/null || true
  cp $raw-counts.txt $baseDir/ena-counts/
  """
}

workflow ena {
  emit: data
  main:
    Channel.fromPath('files/import-data/ena/tpa-urls.txt') | set { urls }
    fetch_metadata(urls) | set { metadata }

    Channel.fromList([
      ['wgs', "$params.databases.ena.remote/wgs/"],
      ['tls', "$params.databases.ena.remote/tls/"],
      ['tsa', "$params.databases.ena.remote/tsa/"],
    ]) \
    | list_subdirs \
    | flatMap { name, listing ->
        // Greedy bin-pack subdirs into batches by inode-size weight (a cheap
        // proxy for entry count), so a few huge subdirs don't make some tasks
        // far heavier than others. A subdir larger than the budget lands in a
        // batch of its own.
        long budget = params.databases.ena.subdir_batch_bytes as long
        def batches = []
        def cur = []
        long curBytes = 0
        listing.readLines().findAll { it.trim() }.each { line ->
          def (sz, dir) = line.split('\t', 2)
          long bytes = sz as long
          if (cur && curBytes + bytes > budget) {
            batches << [name, cur]
            cur = []
            curBytes = 0
          }
          cur << dir
          curBytes += bytes
        }
        if (cur) batches << [name, cur]
        return batches
      } \
    | set { subdir_batches }

    Channel.fromList([
      ['con', ["$params.databases.ena.remote/con/"]],
      ['std', ["$params.databases.ena.remote/std/"]],
    ]) \
    | mix( subdir_batches ) \
    | fetch_directory \
    | flatten \
    | combine(metadata) \
    | process_file \
    | set { data }
}
