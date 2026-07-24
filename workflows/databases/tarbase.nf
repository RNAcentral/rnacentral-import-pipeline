workflow tarbase {
  main:
  remotes = channel.fromList( params.databases.tarbase.remotes )
  remotes | fetch | parse | set { data }

  emit: data
}

process fetch {
  errorStrategy 'retry'
  maxRetries 10

  input:
  val remote
  output:
  path('*.tsv')

  when: { params.databases.tarbase?.run }

  script:
  """
  wget -O tarbase.tsv.gz ${remote}
  gzip -d tarbase.tsv.gz
  """
}

process parse {
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  input:
  path tsv_file

  output:
  path('*.csv')

  when: { params.databases.tarbase?.run }

  script:
  """
  rnac tarbase parse ${tsv_file} .
  """
}
