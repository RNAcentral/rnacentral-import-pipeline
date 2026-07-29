process pombase {
  output:
  path('*.csv')

  when: params.databases.pombase?.run

  script:
  """
  wget -O pombase.json ${params.databases.pombase.remote}
  rnac pombase parse pombase.json .
  """
}
