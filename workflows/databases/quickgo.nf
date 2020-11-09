process fetch_and_process {
  when { params.databases.quickgo.run }

  output:
  path('*.csv')

  """
  scp $params.databases.quickgo.remote data.gpa.gz
  gzip -d data.gpa
  rnac external quickgo data.gpa .
  """
}

workflow quickgo {
  emit: data
  main:
    fetch_and_process | set { data }
}
