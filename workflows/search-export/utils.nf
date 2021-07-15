process query {
  tag { "${query.baseName}" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(max_count)
  path(query)

  output:
  path("${query.baseName}.json")

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > raw.json
  search-export group ${query.baseName} raw.json $max_count ${query.baseName}.json
  """
}
