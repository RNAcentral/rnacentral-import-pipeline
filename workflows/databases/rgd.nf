process rgd {
  when { params.databases.rgd.run }

  output:
  path('*.csv')

  script:
  """
  wget -O sequences.fasta.gz $params.databases.rgd.sequences
  wget -O genes.txt $params.databases.rgd.genes

  rnac rgd parse sequences.fasta.gz genes.txt .
  """
}
