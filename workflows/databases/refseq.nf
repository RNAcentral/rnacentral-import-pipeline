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
  input:
  path(data)

  output:
  path('.csv')

  """
  rnac external refseq $data .
  """
}


workflow refseq {
  emit: data
  main:
    fetch | flatten | parse | set { data }
}
