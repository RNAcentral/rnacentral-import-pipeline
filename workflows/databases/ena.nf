process fetch_single_files {
  when { params.databases.ena.run }

  output:
  path("**/*.ncr.gz")

  script:
  def remote = params.databases.ena.remote
  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$remote/con-std_latest/" con-std

  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$remote/tls/public/" tls

  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$remote/tsa/public/" tsa

  find . -type f -empty -delete
  find . -type f | xargs -I {} gzip --quiet  -l {} | awk '{ if (\$2 == 0) print \$4 }' | xargs -I {} rm {}.gz
  """
}

process find_wgs_directories {
  memory '2GB'

  when { params.databases.ena.run }

  output:
  path('paths')

  script:
  def remote = params.databases.ena.remote
  """
  find $remote/wgs/public -mindepth 1 -maxdepth 1 > paths
  """
}

process fetch_wgs_directories {
  input:
  path(to_fetch)

  output:
  path('files/*')

  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$to_fetch" raw

  mkdir files
  mkdir chunks
  find raw -type f > ids
  split -n l/4000 ids chunks/chunk-
  find chunks/ -type f | xargs -I {} group-wgs {} files
  """
}

process fetch_metadata {
  when { params.databases.ena.run }

  input:
  path(urls)

  output:
  path('tpa.tsv'), emit: tpa
  path('model-lengths.csv'), emit: model_lengths

  """
  cat $urls | xargs -I {} wget -O - {} >> tpa.tsv
  cmstat \$RIBODIR/models/ribo.0p20.extra.cm | grep -v '^#' | awk '{ printf("%s,%d\n", \$2, \$6); }' > model_lengths.csv
  """
}

process process_file {
  tag { "$raw" }

  input:
  tuple path(raw), path(tpa), path(model_lengths)

  output:
  path('*.csv')

  """
  zcat $raw > sequences.dat
  ena2fasta.py sequences.dat sequences.fasta
  /rna/ribovore/ribotyper.pl sequences.fasta ribotyper-results
  rnac ena parse sequences.dat $tpa ribotyper-results $model_lengths .
  """
}

workflow ena {
  emit: data
  main:
    Channel.fromPath('files/import-data/ena/tpa-urls.txt') \
    | fetch_metadata \
    | set { meta }

    find_wgs_directories \
    | splitCsv \
    | map { row -> row[0] } \
    | fetch_wgs_directories \
    | set { wgs_files }

    fetch_single_files \
    | mix(wgs_files) \
    | flatten \
    | combine(meta.out.tpa, meta.out.model_lengths) \
    | process_file \
    | set { data }
}
