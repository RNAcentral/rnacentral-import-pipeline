process fetch {
  queue 'datamover'
  container ''

  output:
  path('circatlas.json')

  script:
  """
  cp ${params.databases.circatlas.remote} circatlas.json
  """
}

process parse {
  memory { params.databases.circatlas.process.directives.memory }

  input:
  path(data)

  output:
  path('*.csv')

  script:
  """
  rnac circatlas parse $data .
  """
}

workflow circatlas {
  emit: data
  main:
  fetch | parse | set { data }
}
