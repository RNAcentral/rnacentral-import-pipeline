process lncbook {
  when: { params.databases.lncbook.run }

  output:
  path('*.csv')

  """
  wget -O lncbook.json ${params.databases.lncbook.remote}
  rnac external lncbook lncbook.json .
  """
}
