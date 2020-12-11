process fetch {
  when { params.databases.silva.run }

  output:
  path('*.rnac')

  """
  wget $params.databases.silva.remote
  gzip -d '*.gz'
  """
}

process parse {
  tag { "$raw.name" }

  input:
  path(raw)

  output:
  path('*.csv')

  """
  rnac external silva $raw .
  """
}

workflow silva {
  emit: data
  main:
    fetch | flatten | parse | set { data }
}
