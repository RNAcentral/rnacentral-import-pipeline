process fetch_data {
  when { params.databases.gtrnadb.run }
  memory '4GB'

  output:
  path('*.json')

  """
  wget --no-check-certificate $params.databases.gtrnadb.remote
  tar xvf *.tar.gz
  """
}

process process_data {
  tag { "$raw.name" }
  memory '4GB'

  input:
  tuple path(raw), path(tax_info)

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
