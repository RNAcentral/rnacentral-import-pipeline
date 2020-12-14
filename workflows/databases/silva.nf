process fetch {
  when { params.databases.silva.run }

  output:
  path('*.rnac')

  """
  wget $params.databases.silva.remote
  gzip -d *.gz
  """
}

process taxonomy_metadata {
  memory '4GB'
  when { params.databases.silva.run }

  output:
  path('taxonomy.db')

  """
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  tar xvf new_taxdump.tar.gz
  mkdir taxdump
  mv *.dmp taxdump
  rnac silva index-taxonomy taxdump taxonomy.db
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
  emit: data
  main:
    fetch \
    | flatten \
    | combine(taxonomy_metadata())
    | parse \
    | set { data }
}
