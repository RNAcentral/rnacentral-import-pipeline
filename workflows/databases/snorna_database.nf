process snorna_database {
  when: { params.databases.snorna_database.run }

  output:
  path('*.csv')

  """
  scp ${params.databases.snorna_database.remote} snorna_database.json 
  rnac snorna_database parse snorna_database.json .
  """
}
