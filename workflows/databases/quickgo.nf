process quickgo {
  when { params.databases.quickgo.run }
  memory { params.databases.quickgo.memory }

  output:
  path('*.csv')

  """
  scp $params.databases.quickgo.remote data.gpa.gz
  gzip -d data.gpa
  rnac quickgo parse data.gpa .
  """
}
