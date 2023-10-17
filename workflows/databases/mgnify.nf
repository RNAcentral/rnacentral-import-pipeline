process mgnify_fetch {
  queue 'datamover'
  container ''
  output:
    path('*.json')

  """
  cp /nfs/ftp/public/databases/metagenomics/rnacentral/mgnify_genomes/* .
  """
}

process mgnify_parse {
  input:
    path(json)
  output:
    path('*.csv')

  """
  rnac mgnify parse $json .
  """
}


workflow mgnify {
  main:
    if ( params.databases.mgnify.run ) {
      mgnify_fetch | mgnify_parse | set { data }
    }
    else {
      Channel.empty() | set { data }
    }
}
