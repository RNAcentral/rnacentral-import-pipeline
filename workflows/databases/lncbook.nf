process fetch_and_process {
  when: { params.databases.lncbook.run }

  output:
  path('*.csv')

  """
  wget -O lncbook.json ${params.databases.lncbook.remote}
  rnac external lncbook lncbook.json .
  """
}

workflow lncbook {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}


