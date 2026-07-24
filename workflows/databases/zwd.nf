process zwd {
  input:
  path(context)

  output:
  path('*.csv')

  when: { params.databases.zwd?.run }

  script:
  """
  wget -O zwd.json $params.databases.zwd.remote
  rnac zwd parse $context zwd.json .
  """
}
