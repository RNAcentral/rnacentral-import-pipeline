process five_s_rrnadb {
  output:
  path('*.csv')

  when: { params.databases["5srrnadb"]?.run }

  script:
  """
  scp ${params.databases["5srrnadb"].remote} 5s.json
  rnac 5srrnadb parse 5s.json .
  """
}
