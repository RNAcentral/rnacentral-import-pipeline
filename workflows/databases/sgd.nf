process fetch_and_process {
  when: { params.databases.sgd.run }

  output:
  path('*.csv')

  """
  wget -O sgd.json ${params.databases.sgd.remote}
  rnac external sgd sgd.json .
  """
}

workflow sgd {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}

