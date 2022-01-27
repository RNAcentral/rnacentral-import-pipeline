process ribovision {
  output:
  path("*.csv")

  """
  curl ${params.databases.ribovision.remote} > raw.html
  rnac ribovision parse raw.html
  """
}
