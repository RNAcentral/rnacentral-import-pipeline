process fetch {
  when { params.databases.silva.run }

  output:
  path('*.rnac.gz')

  """
  wget $params.silva.remote
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
    fetch | parse | set { data }
}
