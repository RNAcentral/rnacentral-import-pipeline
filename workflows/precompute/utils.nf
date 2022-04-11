process query {
  tag { "${query.baseName}" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(max_count)
  path(query)

  output:
  path("${query.baseName}.json")

  """
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=$params.precompute.tablename \
    -f $query \
    "$PGDATABASE" > raw.json
  precompute metadata group ${query.baseName} raw.json $max_count ${query.baseName}.json
  """
}
