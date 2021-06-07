process fetch_directory {
  tag { "$name" }
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

  input:
  tuple val(name), val(remote)

  output:
  path "$name/*.ncr.gz", emit: single_files
  path "$name/*.gz", emit: tar_files

  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --include'**/*.tar'
    --exclude='*.fasta.gz' \
    "$remote" "$name"

  find $name -type f -empty -delete
  find $name -type f -name '*.gz' | xargs -I {} gzip --quiet -l {} | awk '{ if (\$2 == 0) print \$4 }' | xargs -I {} rm {}.gz
  """
}

process expand_tar_files {
  tag { "${to_fetch.name}" }

  input:
  path(tar_file)

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

    Channel.fromList([
      ['con', "$params.databases.ena.remote/con/"],
      ['std', "$params.databases.ena.remote/std/"],
      ['tls', "$params.databases.ena.remote/tls/"],
      ['tsa', "$params.databases.ena.remote/tsa/"],
      ['wgs', "$params.databases.ena.remote/wgs/"],
    ]) \
    | fetch_directory

    fetch_directory.tar_files | expand_tar_files | set { tar_sequences }

    fetch_directory.single_files \
    | mix(tar_sequences) \
    | flatten \
    | combine(metadata) \
    | process_file \
    | set { data }
}
