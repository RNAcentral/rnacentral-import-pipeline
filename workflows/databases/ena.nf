process fetch_directory {
  tag { "$name" }
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

  input:
  tuple val(name), val(remote)

  output:
  path("${name}-chunks/*.ncr")

  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --include='**/*.tar' \
    --exclude='*.fasta.gz' \
    "$remote" "copied"

  find copied -type f -empty -delete
  find copied -type f -name '*.gz' | xargs -I {} gzip --quiet -l {} | awk '{ if (\$2 == 0) print \$4 }' | xargs -I {} rm {}.gz

  mkdir $name
  pushd $name
  find . -type f -name '*.tar' | xargs -I {} tar -xvf {}
  popd
  find copied -type f -name '*.ncr.gz' | xargs zcat > ${name}.ncr

  mkdir $name-chunks
  split-ena --max-sequences ${params.databases.ena.max_sequences} ${name}.ncr ${name}-chunks
  """
}

process fetch_metadata {
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

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
      ['con', "$params.databases.ena.remote/con/"],
      ['std', "$params.databases.ena.remote/std/"],
      ['tls', "$params.databases.ena.remote/tls/"],
      ['tsa', "$params.databases.ena.remote/tsa/"],
      ['wgs', "$params.databases.ena.remote/wgs/"],
    ]) \
    | fetch_directory \
    | flatten \
    | combine(metadata) \
    | process_file \
    | set { data }
}
