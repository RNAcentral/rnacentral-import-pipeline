process fetch_data {
  when: { params.databases.mgnify.run }
  queue 'datamover'
  output:
    path('*.json')

  """
  cp /nfs/ftp/public/databases/metagenomics/rnacentral/mgnify_genomes/* .
  """
}

process parse {
  when: { params.databases.mgnify.run }
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
