process fetch {
  when { params.databases.modomics.run }
  queue 'datamover'
  container ''

  output:
  path('modomics.json')

  """
  cp ${params.databases.modomics.remote} modomics.json
  """
}

process parse {
  when { params.databases.modomics.run }

  input:
  path(data)

  output:
  path('*.csv')

  """
  rnac modomics parse $data .
  """
}

workflow modomics {
  emit: data
  main:
  fetch | parse | set { data }
}
