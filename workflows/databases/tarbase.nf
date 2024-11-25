workflow tarbase {
  emit: data
  remotes = channel.of( params.databases.tarbase.remotes )
  remotes | fetch | parse | set { data }
}

process fetch {
  when: { params.databases.tarbase.run }

  input:
  val remote
  output:
  path('*.tsv')

  """
  wget -O tarbase.tsv.gz ${remote}
  gzip -d tarbase.tsv.gz
  """
}

process parse {
  when: { params.databases.tarbase.run }

  input:
  path tsv_file

  output:
  path('*.csv')

  """
  rnac tarbase parse-tsv ${tsv_file} .
  """
}
