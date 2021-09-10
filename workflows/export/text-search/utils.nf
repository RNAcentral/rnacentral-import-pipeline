process fetch {
  tag { "${query.baseName}" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(max_count)
  path(query)

  output:
  tuple val(query.baseName), path("raw.json"), val(max_count)

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > raw.json
  """
}

process group {
  tag { "${query.baseName}" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(name), val("raw.json"), val(max_count)

  output:
  path("${name}.json")

  """
  search-export group ${name} raw ${max_count} ${name}.json
  """
}

workflow query {
  take:
    max_count
    query
  emit: data
  main:
    fetch(max_count, query) | group | set { data }
}
