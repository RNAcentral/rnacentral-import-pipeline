process fetch_data {
  when { params.databases.gtrnadb.run }

  output:
  path('*.json')

  """
  scp $params.databases.gtrnadb.remote .
  gzip -d *.json.gz
  gzip -d bacteria_tRNAs.tar.gz
  mv bacteria_tRNAs.tar bacteria_tRNAs.json
  ls *.tar.gz | xargs -I {} --verbose tar -xzvf {}
  """
}

process process_data {
  input:
  path('raw.json')

  output:
  path('*.csv')

  """
  rnac external gtrnadb raw.json .
  """
}

workflow gtrnadb {
  emit: data
  main:
    fetch_data | process_data | set { data }
}
