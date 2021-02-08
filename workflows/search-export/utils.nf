process query {
  tag { "${query.baseName}" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  val(_flag)
  path(query)

  output:
  path("${query.baseName}.json")

  """
  psql \
    --variable ON_ERROR_STOP=1 \
    -f $query \
    "$PGDATABASE" > ${query.baseName}.json
  """
}
