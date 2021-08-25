process process_data {
  output:
  path("*.csv")

  """
  cp $params.databases.psicquic.remote raw.tsv
  rnac psicquic parse raw.tsv .
  """
}

workflow psicquic {
  emit: data
  main:
    process_data
}

workflow {
  psicquic()
}
