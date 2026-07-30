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
  tuple path(raw), path(tpa), path(model_lengths)

  output:
  path('*.{csv,parquet}')

  script:
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
    | combine(metadata) \
    | process_file \
    | set { data }

  emit: data
}
