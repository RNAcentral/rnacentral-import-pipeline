process fetch_and_process {
  when: { params.databases.flybase.run }

  output:
  path('*.csv')

  """
  wget -O - ${params.databases.flybase.remote} | gzip -d > flybase.json
  rnac external flybase flybase.json .
  """
}

workflow flybase {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
