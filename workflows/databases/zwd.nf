process zwd {
  when: { params.databases.zwd.run }

  input:
  path(context)

  output:
  path('*.{csv,parquet}')

  """
  wget -O zwd.json $params.databases.zwd.remote
  rnac zwd parse $context zwd.json .
  """
}
