process fetch_and_process {
  when { params.databases.mirbase.run }

  output:
  path('*.csv')

  """
  scp $params.databases.mirbase.remote mirbase.json
  rnac external mirbase mirbase.json .
  """
}

workflow mirbase {
  emit: data
  main:
    fetch_and_process | set { data }
}
