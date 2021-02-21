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
  precompute group ${query.baseName} raw.json $max_count ${query.baseName}.json
  """
}

process fetch_release_info {
  when { params.precompute.run }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory '10GB'

  input:
  val(_flag)
  path(query)

  output:
  path('data.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > raw
  sort -t , -nk1,1 -u raw > sorted.csv
  precompute max-release sorted.csv data.csv
  """
}
