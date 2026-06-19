process psicquic {
  when: { params.databases.psicquic?.run }

  output:
  path("*.csv")

  script:
  """
  cp $params.databases.psicquic.remote raw.tsv
  rnac psicquic parse raw.tsv .
  """
}
