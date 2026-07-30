process ribovision {
  output:
  path('*.{csv,parquet}')

  when: params.databases.ribovision?.run

  script:
  """
  curl ${params.databases.ribovision.remote} > raw.html
  rnac ribovision parse raw.html
  """
}
