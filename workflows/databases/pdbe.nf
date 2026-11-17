process pdbe_fetch {
  queue 'datamover'
  container ''

  output:
  path('pdb_full_region.txt.gz')

  when: params.databases.pdb?.run

  script:
  """
  cp /nfs/ftp/public/databases/Rfam/.preview/pdb_full_region.txt.gz .
  """
}

process pdbe_generate {
  errorStrategy 'retry'
  maxRetries 5

  input:
  path(regions)

  output:
  path('*.{csv,parquet}')

  when: params.databases.pdb?.run

  script:
  """
  gzip -cd $regions | awk 'BEGIN {OFS = FS = "\t" } \$11 == 1 { print \$2, \$3} ' | sort -u  > rfam_hit_ids
  rnac pdb generate --override-chains=rfam_hit_ids .
  """
}

workflow pdbe {
  main:
  pdbe_fetch | pdbe_generate | set { data }

  emit: data
}
