process fetch_and_process {
  when: { params.databases.zwd.run }

  output:
  path('*.csv')

  """
  curl $params.databases.zwd.remote > zwd.json
  rnac external zwd zwd.json .
  """
}

workflow zwd {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
