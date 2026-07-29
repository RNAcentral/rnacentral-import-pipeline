process fetch {
  output:
  path('*.gbff')

  when: params.databases.refseq?.run

  script:
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

  script:
  """
  rnac refseq parse $data .
  """
}


workflow refseq {
  main:
    fetch | flatten | parse | set { data }
  emit: data
}
