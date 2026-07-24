process mgnify_fetch {
  queue 'datamover'
  container ''
  output:
    path('*.json')

  when: { params.databases.mgnify?.run }

  script:
  """
  cp /nfs/ftp/public/databases/metagenomics/rnacentral/mgnify_genomes/* .
  """
}

process mgnify_parse {
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  input:
    path(json)
  output:
    path('*.csv')

  script:
  """
  rnac mgnify parse $json .
  """
}


workflow mgnify {
  main:
    if ( params.databases.mgnify?.run ) {
      mgnify_fetch | flatten | mgnify_parse | set { data }
    }
    else {
      Channel.empty() | set { data }
    }

  emit: data
}
