process fetch_and_process {
  when: { params.databases.tarbase.run }

  output:
  path('*.csv')

  """
  cp ${params.databases.tarbase.remote} tarbase.json
  rnac external tarbase tarbase.json .
  """
}

workflow tarbase {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}

