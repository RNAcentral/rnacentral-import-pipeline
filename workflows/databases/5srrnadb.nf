process five_s_rrnadb {
  when: { params.databases["5srrnadb"].run }

  output:
  path('*.csv')

  """
  scp ${params.databases["5srrnadb"].remote} 5s.json
  rnac external 5srrnadb 5s.json .
  """
}
