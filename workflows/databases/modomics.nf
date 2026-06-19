process fetch {
  queue 'datamover'
  container ''

  output:
  path('modomics.json')

  when: { params.databases.modomics.run }

  script:
  """
  cp ${params.databases.modomics.remote} modomics.json
  """
}

process parse {
  input:
  path(data)

  output:
  path('*.csv')

  script:
  """
  rnac modomics parse $data .
  """
}

workflow modomics {
  emit: data
  main:
  fetch | parse | set { data }
}
