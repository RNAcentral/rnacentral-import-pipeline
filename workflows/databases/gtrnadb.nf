process fetch_data {
  when { params.databases.gtrnadb.run }

  output:
  path('*.json')

  """
  rnac gtrnadb urls-for $params.databases.gtrnadb.remote | xargs -I {} wget {}
  gzip -d *.json.gz
  find . -name '*.tar.gz' | xargs -I {} --verbose tar -xzvf {}
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
    fetch_data | flatten | process_data | set { data }
}
