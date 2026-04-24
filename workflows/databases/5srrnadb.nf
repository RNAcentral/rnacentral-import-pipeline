process five_s_rrnadb {
  when: { params.databases["5srrnadb"].run }

  output:
  path('*.{csv,parquet}')

  """
  scp ${params.databases["5srrnadb"].remote} 5s.json
  rnac 5srrnadb parse 5s.json .
  """
}
