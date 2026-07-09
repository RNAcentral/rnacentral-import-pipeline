process ribovision {
  when: params.databases.ribovision?.run

  output:
  path("*.csv")

  script:
  """
  curl ${params.databases.ribovision.remote} > raw.html
  rnac ribovision parse raw.html
  """
}
