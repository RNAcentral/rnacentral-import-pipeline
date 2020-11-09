process fetch_and_process {
  when: { params.databases.lncipedia.run }

  output:
  path('*.csv')

  """
  curl ${params.databases.lncipedia.remote} > lncipedia.json
  rnac external lncipedia lncipedia.json .
  """
}

workflow lncipedia {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
