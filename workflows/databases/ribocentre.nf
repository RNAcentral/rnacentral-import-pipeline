
process fetch_ribocentre {
  memory 1.GB

  output:
  path("ribocentre.json")

  when: { params.databases.ribocentre?.run }

  script:
  """
  wget -O ribocentre.json ${params.databases.ribocentre.remote}
  """
}

process parse_ribocentre {
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  input:
  path ribocentre_json
  output:
  path("*.csv")

  when: { params.databases.ribocentre?.run }

  script:
  """
  rnac ribocentre parse $ribocentre_json
  """
}


workflow ribocentre {
  main:
  fetch_ribocentre | parse_ribocentre | set { data }

  emit: data
}
