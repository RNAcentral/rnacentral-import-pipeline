process fetch {
  queue 'datamover'
  container ''

  output:
  path('circatlas.json')

  when: params.databases.circatlas?.run

  script:
  """
  cp ${params.databases.circatlas.remote} circatlas.json
  """
}

process parse {
  memory '16 GB'

  input:
  path(data)

  output:
  path('*.csv')

  when: params.databases.circatlas?.run

  script:
  """
  rnac circatlas parse $data .
  """
}

workflow circatlas {
  main:
  fetch | parse | set { data }
  emit: data
}
