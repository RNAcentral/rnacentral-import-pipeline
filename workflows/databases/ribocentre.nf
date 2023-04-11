process ribocentre {
  when: params.databases.ribocentre.run

  output:
  path("*.csv")

  """
  rnac ribocentre parse ${params.databases.ribocentre.remote}
  """
}
