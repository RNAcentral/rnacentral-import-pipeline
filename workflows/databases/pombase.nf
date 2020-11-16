process fetch_and_process {
  when: { params.databases.pombase.run }

  output:
  path('*.csv')

  """
  wget -O pombase.json ${params.databases.pombase.remote}
  rnac external pombase pombase.json .
  """
}

workflow pombase {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
