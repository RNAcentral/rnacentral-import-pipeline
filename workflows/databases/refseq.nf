process fetch {
  when { params.databases.refseq.run }

  output:
  path('*.gbff')

  """
  wget "$params.databases.refseq.remote"
  gzip -d *.gbff.gz
  """
}

process parse {
  tag { "$data.name" }

  input:
  path(data)

  output:
  path('*.csv')

  """
  rnac refseq parse $data .
  """
}


workflow refseq {
  emit: data
  main:
    fetch | flatten | parse | set { data }
}
