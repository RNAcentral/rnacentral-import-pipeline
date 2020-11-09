process fetch_and_process {
  when { params.databases.snodb.run }

  output:
  path('*.csv')

  """
  scp $params.databases.snodb.remote snodb.json
  rnac external snodb snodb.json .
  """
}

workflow snodb {
  emit: data
  main:
    fetch_and_process | set { data }
}
