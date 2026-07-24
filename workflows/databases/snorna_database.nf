process snorna_database {
  output:
  path('*.csv')

  when: { params.databases.snorna_database?.run }

  script:
  """
  scp ${params.databases.snorna_database.remote} snorna_database.json 
  rnac snorna_database parse snorna_database.json .
  """
}
