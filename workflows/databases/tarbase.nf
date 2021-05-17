process tarbase {
  when: { params.databases.tarbase.run }

  output:
  path('*.csv')

  """
  cp ${params.databases.tarbase.remote} tarbase.json
  rnac tarbase parse tarbase.json .
  """
}
