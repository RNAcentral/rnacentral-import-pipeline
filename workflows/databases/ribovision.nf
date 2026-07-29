process ribovision {
  output:
  path("*.csv")

  when: params.databases.ribovision?.run

  script:
  """
  curl ${params.databases.ribovision.remote} > raw.html
  rnac ribovision parse raw.html
  """
}
