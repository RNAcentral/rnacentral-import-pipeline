process fetch_single_files {
  tag { "$name" }
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

  input:
  tuple val(name), val(remote)

  output:
  path("**/*.ncr.gz")

  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$remote" "$name"

  find $name -type f -empty -delete
  find $name -type f | xargs -I {} gzip --quiet -l {} | awk '{ if (\$2 == 0) print \$4 }' | xargs -I {} rm {}.gz

  find $name -type f |\
  xargs -I {} zgrep -Hc '^ID' {} |\
  awk -F ':' '{ if (\$2 > 100000) print \$1 }' |\
  xargs -I {} gzip -d {}

  mkdir large-splits
  find $name -type f -name '*.ncr' | xargs -I {} split-ena --max-sequences $params.databases.ena.large_chunk_size --remove-file {} large-splits
  find large-splits -type f | xargs -I {} gzip {}
  """
}

process find_wgs_directories {
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

  output:
  path('paths')

  script:
  def remote = params.databases.ena.remote
  """
  find $remote/wgs/public -mindepth 1 -maxdepth 1 -not -empty > paths
  """
}

process fetch_wgs_directories {
  tag { "${file(to_fetch).name}" }
  clusterOptions '-sp 90'

  input:
  val(to_fetch)

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
  split -n l/60 ids chunks/chunk-
  find chunks -type f -empty -delete
  find chunks/ -type f | xargs -I {} group-wgs {} "${file(to_fetch).name}" files
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
  memory '4GB'
  tag { "$raw" }

  input:
  tuple path(raw), path(tpa), path(model_lengths)

  output:
  path('*.csv')

  """
  zcat $raw > sequences.dat
  ena2fasta.py sequences.dat sequences.fasta
  if [[ -e sequences.fasta ]]; then
    /rna/ribovore/ribotyper.pl sequences.fasta ribotyper-results
  else
    mkdir ribotyper-results
  fi
  rnac ena parse --counts $raw-counts.txt sequences.dat $tpa ribotyper-results $model_lengths .

  mkdir $baseDir/ena-counts 2>/dev/null || true
  cp $raw-counts.txt $baseDir/ena-counts/
  """
}

workflow ena {
  emit: data
  main:
    Channel.fromPath('files/import-data/ena/tpa-urls.txt') | set { urls }

    find_wgs_directories \
    | splitCsv \
    | map { row -> row[0] } \
    | fetch_wgs_directories \
    | set { wgs_files }

    Channel.fromList([
      ['con-std', "$params.databases.ena.remote/con-std_latest/"],
      ['tls', "$params.databases.ena.remote/tls/public/"],
      ['tsa', "$params.databases.ena.remote/tsa/public/"],
    ]) \
    | fetch_single_files \
    | mix(wgs_files) \
    | flatten \
    | combine(fetch_metadata(urls)) \
    | process_file \
    | set { data }
}
