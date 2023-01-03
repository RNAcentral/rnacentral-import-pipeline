process ribocentre {
  when: params.databases.ribocentre.run

  output:
  path("*.csv")

  """
  rnac ribovision parse ${params.databases.ribocentre.remote}
  """
