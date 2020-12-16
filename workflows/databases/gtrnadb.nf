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
  tag { "$raw.name" }

  input:
  tuple path(tax_info), path(raw)

  output:
  path('*.csv')

  """
  rnac gtrnadb parse $tax_info $raw .
  """
}

workflow gtrnadb {
  take: tax_info
  emit: data
  main:
    fetch_data \
    | flatten \
    | combine(tax_info) \
    | process_data \
    | set { data }
}
