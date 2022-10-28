process taxonomy {
  memory '2GB'
  errorStrategy 'retry'

  output:
  path('*.csv')

  """
  wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5
  md5sum -c new_taxdump.tar.gz.md5
  tar xvf new_taxdump.tar.gz
  mkdir taxdump
  mv *.dmp taxdump
  rnac ncbi taxonomy taxdump
  """
}
