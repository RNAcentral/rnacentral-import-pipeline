
process fetch_ribocentre {
  memory 1.GB

  when: { params.databases.ribocentre.run }
  input:
  val ribocentre_remote
  output:
  path("ribocentre.json")

  """
  wget -O ribocentre.json ${params.databases.ribocentre.remote}
  """
}

process parse_ribocentre {
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  when: { params.databases.ribocentre.run }
  input:
  path ribocentre_json
  output:
  path("*.csv")

  """
  rnac ribocentre parse $ribocentre_json
  """
}


workflow ribocentre {
  emit: data

  main:
  fetch_ribocentre | parse_ribocentre | set { data }
}
