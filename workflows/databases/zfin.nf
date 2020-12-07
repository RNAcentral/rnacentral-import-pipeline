process zfin {
  when: { params.databases.zfin.run }

  output:
  path('*.csv')

  """
  rnac fetch zfin$params.databases.zfin.remote zwd.json
  rnac external zfin zfin.json .
  """
}
