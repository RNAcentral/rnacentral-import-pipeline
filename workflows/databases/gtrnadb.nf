process fetch_data {
  queue 'datamover'
  container ''

  output:
  path('*.json')

  when: params.databases.gtrnadb?.run

  script:
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
  path('*.{csv,parquet}')

  script:
  """
  rnac gtrnadb parse $tax_info $raw .
  """
}

workflow gtrnadb {
  take: tax_info
  main:
    fetch_data \
    | flatten \
    | combine(tax_info) \
    | process_data \
    | set { data }
  emit: data
}
