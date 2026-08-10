process zwd {
  input:
  path(context)

  output:
  path('*.{csv,parquet}')

  when: params.databases.zwd?.run

  script:
  """
  wget -O zwd.json $params.databases.zwd.remote
  rnac zwd parse $context zwd.json .
  """
}
