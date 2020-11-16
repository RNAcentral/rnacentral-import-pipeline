process fetch_data {
  when { params.databases.ena.run }

  output:
  path("*.ncr"), emit: single_file
  path("wgs/public/*"), emit: directories

  script:
  def remote = params.databases.ena.remote
  """
  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='*/' \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$remote/con-std_latest" "$remote/tls/public" "$remote/tsa/public" .

  rsync \
    -avPL \
    --prune-empty-dirs \
    --include='**/*.ncr.gz' \
    --exclude='*.fasta.gz' \
    "$remote/wgs/public/*" .
  """
}

process fetch_tpa {
  when { params.databases.ena.run }

  input:
  path(urls)

  output:
  path('tpa.tsv')

  """
  cat $urls | xargs -I {} wget -O - {} >> tpa.tsv
  """
}

process process_data {
  tag { "$raw" }

  input:
  tuple path(raw), path(tpa)

  output:
  path('*.csv')

  """
  ena2fasta.py $raw sequences.fasta
  ribotyper.pl sequences.fasta ribotyper-results
  rnac external ena $raw $tpa ribotyper-results
  """
}

workflow ena {
  emit: data
  main:
    Channel.fromPath('files/import-data/ena/tpa-urls.txt') | fetch_tpa | set { tpa }

    fetch_data()

    mix(
      fetch_data.out.single_file,
      fetch_data.out.directories,
    ) \
    | combine(tpa) \
    | process_data \
    | set { data }
}
