process rgd {
  when { params.databases.rgd.run }

  output:
  path('*.csv')

  """
  wget -O sequences.fasta.gz $params.databases.rgd.sequences
  wget -O genes.txt $params.databases.rgd.genes
  gzip -d sequences.fasta.gz

  rnac external rgd sequences.fasta genes.txt .
  """
}
