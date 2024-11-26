workflow tarbase {
  emit: data
  remotes = channel.fromList( params.databases.tarbase.remotes )
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
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  when: { params.databases.tarbase.run }

  input:
  path tsv_file

  output:
  path('*.csv')

  """
  rnac tarbase parse ${tsv_file} .
  """
}
