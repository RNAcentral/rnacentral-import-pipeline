process ribocentre {
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  when: { params.databases.ribocentre.run }
  output:
  path("*.csv")

  """
  rnac ribocentre parse ${params.databases.ribocentre.remote}
  """
}
