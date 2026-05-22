process pombase {
  when: { params.databases.pombase.run }

  output:
  path('*.csv')

  script:
  """
  wget -O pombase.json ${params.databases.pombase.remote}
  rnac pombase parse pombase.json .
  """
}
