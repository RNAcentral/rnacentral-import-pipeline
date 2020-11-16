process fetch_and_process {
  when: { params.databases.zfin.run }

  output:
  path('*.csv')

  """
  rnac fetch zfin$params.databases.zfin.remote zwd.json
  rnac external zfin zfin.json .
  """
}

workflow zfin {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
