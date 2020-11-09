process fetch_and_process {
  when: { params.databases.snorna_database.run }

  output:
  path('*.csv')

  """
  scp ${params.databases.snorna_database.remote} snorna_database.json 
  rnac external snorna_database snorna_database.json .
  """
}

workflow snorna_database {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
