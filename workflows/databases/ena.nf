process fetch_directory {
  tag { "$name" }
  when { params.databases.ena.run }
  clusterOptions '-sp 100'

  input:
  tuple val(name), val(remote)

  output:
  path "$name/*.{ncr.gz,tar}"

  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --include='**/*.tar' \
    --exclude='*.fasta.gz' \
    "$remote" "$name"

  find $name -type f -empty -delete
  find $name -type f -name '*.gz' | xargs -I {} gzip --quiet -l {} | awk '{ if (\$2 == 0) print \$4 }' | xargs -I {} rm {}.gz
  """
}

process expand_tar_files {
  tag { "${to_fetch.name}" }
  clusterOptions '-sp 90'

  input:
  path(tar_file)

  output:
  path("${tar_file.simpleName}/*.ncr.gz")

  """
  tar -xvf "$tar_file"
  mv ${tar_file.simpleName} expanded
  find "expanded" -name '*.ncr.gz' | xargs zcat > all-sequences.ncr
  mkdir ${tar_file.simpleName}
  split-ena --max-sequences ${params.databases.ena.max_sequences} all-sequences.ncr ${tar_file.simpleName}
  find ${tar_file.simpleName} -type f | xargs -I {} gzip {}
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
    | fetch_directory \
    | flatten \
    | branch {
      tar: it.getExtension() == 'tar'
      gz: it.getExtension() == 'gz'
      other: true
    } \
    | set { files }

    files.tar | expand_tar_files | set { tar_sequences }
    files.other.map { error("Unknown fetched file ${it}") }

    files.gz \
    | mix(tar_sequences) \
    | flatten \
    | combine(metadata) \
    | process_file \
    | set { data }
}
