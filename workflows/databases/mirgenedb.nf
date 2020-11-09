process fetch_and_process {
  when { params.databases.mirgenedb.run }

  output:
  path('*.csv')

  """
  scp $params.databases.mirgenedb.remote mirgenedb.json
  rnac external mirgenedb mirgenedb.json .
  """
}

workflow mirgenedb {
  emit: data
  main:
    fetch_and_process | set { data }
}
