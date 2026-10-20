process quickgo_get {
  queue 'datamover'
  container ''

  output:
  path('data.gpa')

  script:
  """
  scp $params.databases.quickgo.remote data.gpa.gz
  gzip -d data.gpa.gz
  """
}



process quickgo_parse {
  memory { params.databases.quickgo.memory * task.attempt }
  errorStrategy { task.attempt <= 3 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
  path(data)

  output:
  path('*.{csv,parquet}')

  script:
  """
  rnac quickgo parse $data .
  """
}



workflow quickgo {

  main:
    if ( params.databases.quickgo?.run ) {
      quickgo_get | quickgo_parse | set { data }
    }
    else {
      channel.empty() | set { data }
    }

  emit: data

}
