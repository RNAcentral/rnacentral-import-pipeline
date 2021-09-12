process fetch_data {
  when { params.databases.gtrnadb.run }

  output:
  path('*.json')

  """
  wget $params.databases.gtrnadb.remote
  tar xvf *.tar.gz
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
