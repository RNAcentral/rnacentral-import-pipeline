process fetch_single_files {
  tag { "$name" }
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

  input:
  tuple val(name), val(remote)

  output:
  path("$name/*.ncr.gz")

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
  cp $remote/wgs/*.tar .
  """
}

process expand_wgs_directories {
  tag { "${file(to_fetch).name}" }

  input:
  val(tar_file)

  output:
  path("${tar_file.name}.ncr.gz")

  """
  tar -xvf "$tar_file"
  find . -name '*.ncr.gz' | zcat > "${tar_file.name}.ncr.gz"
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
    fetch_metadata(urls) | set { metadata }

    fetch_wgs_directories \
    | expand_wgs_directories \
    | set { wgs_files }

    Channel.fromList([
      ['con', "$params.databases.ena.remote/con/"],
      ['std', "$params.databases.ena.remote/std/"],
      ['tls', "$params.databases.ena.remote/tls/"],
      ['tsa', "$params.databases.ena.remote/tsa/"],
    ]) \
    | fetch_single_files \
    | mix(wgs_files) \
    | flatten \
    | combine(metadata) \
    | process_file \
    | set { data }
}
