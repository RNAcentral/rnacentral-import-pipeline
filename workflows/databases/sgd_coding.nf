process sgd_coding {
  memory '4GB'

  output:
  path('*.{csv,parquet}')

  when: params.databases.sgd_coding?.run

  script:
  """
  wget -O - ${params.databases.sgd_coding.remote} | gzip -dc > sgd.gff
  rnac sgd parse-coding sgd.gff .
  """
}
