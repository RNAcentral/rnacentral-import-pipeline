process ribocentre {
  when: { params.databases.ribocentre.run }
  memory 8.GB
  output:
  path("*.csv")

  """
  rnac ribocentre parse ${params.databases.ribocentre.remote}
  """
}
