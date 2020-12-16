process fetch {
  when { params.databases.silva.run }

  output:
  path('*.rnac')

  """
  wget $params.databases.silva.remote
  gzip -d *.gz
  """
}

process parse {
  tag { "$raw.name" }

  input:
  tuple path(raw), path(taxonomy)

  output:
  path('*.csv')

  """
  rnac silva parse $raw $taxonomy .
  """
}

workflow silva {
  take: tax_info
  emit: data
  main:
    fetch \
    | flatten \
    | combine(tax_info) \
    | parse \
    | set { data }
}
