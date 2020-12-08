process fetch {
  when { params.databases.silva.run }

  output:
  path('*.rnac.gz')

  """
  wget $params.databases.silva.remote
  """
}


process parse {
  input:
  path(raw)

  output:
  path('*.csv')

  """
  zcat $raw | rnac external silva - .
  """
}


workflow silva {
  emit: data
  main:
    fetch | flatten | parse | set { data }
}
