process quickgo {
  when { params.databases.quickgo.run }

  output:
  path('*.csv')

  """
  scp $params.databases.quickgo.remote data.gpa.gz
  gzip -d data.gpa
  rnac external quickgo data.gpa .
  """
}
