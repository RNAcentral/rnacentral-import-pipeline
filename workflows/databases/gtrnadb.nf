process fetch_data {
  when { params.databases.gtrnadb.run }
  queue 'datamover'
  container ''

  output:
  path('*.json')

  """
  cp $params.databases.gtrnadb.remote/*.json .
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
