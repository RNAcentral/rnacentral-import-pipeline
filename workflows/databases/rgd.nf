process rgd {
  output:
  path('*.{csv,parquet}')

  when: { params.databases.rgd?.run }

  script:
  """
  wget -O sequences.fasta.gz $params.databases.rgd.sequences
  wget -O genes.txt $params.databases.rgd.genes

  rnac rgd parse sequences.fasta.gz genes.txt .
  """
}
