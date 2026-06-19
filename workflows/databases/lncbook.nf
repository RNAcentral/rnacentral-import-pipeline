process lncbook {
  when: { params.databases.lncbook?.run }

  output:
  path('*.csv')

  script:
  """
  wget -O lncbook.json.gz ${params.databases.lncbook.remote}
  gzip -d lncbook.json.gz
  rnac lncbook parse lncbook.json .
  """
}
