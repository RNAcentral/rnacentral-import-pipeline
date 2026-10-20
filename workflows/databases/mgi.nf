process mgi {
  output:
  path('*.{csv,parquet}')

  when: params.databases.mgi?.run

  script:
  """
  wget -O MRK_Sequence.rpt '$params.databases.mgi.remote'
  rnac mgi parse MRK_Sequence.rpt .
  """
}
