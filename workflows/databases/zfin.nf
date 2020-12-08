process zfin {
  when: { params.databases.zfin.run }

  output:
  path('*.csv')

  """
  rnac zfin fetch $params.databases.zfin.remote zwd.json
  rnac zfin parse zfin.json .
  """
}
