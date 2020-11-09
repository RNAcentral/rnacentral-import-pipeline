process fetch_and_process {
  when: { params.databases.lncbase.run }

  output:
  path('*.csv')

  """
  cp ${params.databases.lncbase.remote} lncbase.json
  rnac external lncbase lncbase.json .
  """
}

workflow lncbase {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}

