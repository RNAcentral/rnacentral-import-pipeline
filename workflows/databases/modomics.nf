process fetch {
  queue 'datamover'
  container ''

  output:
  path('modomics.json')

  when: params.databases.modomics?.run

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

  when: params.databases.modomics?.run

  script:
  """
  rnac modomics parse $data .
  """
}

workflow modomics {
  main:
  fetch | parse | set { data }

  emit: data
}
