process pdbe {
  when { params.databases.pdb.run }

  output:
  path('*.csv')

  """
  wget --read-timeout=30 -t 1 -O pdb_full_region.txt.gz http://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz
  gzip -d pdb_full_region.txt.gz
  awk 'BEGIN {OFS = FS = "\t" } \$11 == 1 { print \$2, \$3} ' pdb_full_region.txt | sort -u  > rfam_hit_ids
  rnac pdb generate --override-chains=rfam_hit_ids .
  """
}
