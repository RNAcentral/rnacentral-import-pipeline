params.maxForks = 10

process simple_query {
  maxForks params.maxForks

  input:
  path(query)

  output:
  path('result')

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > result
  """
}

process with_params {
  maxForks params.max_forks

  input:
  tuple path(query), val(parameters)

  output:
  path('raw')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" $parameters "$PGDATABASE" > raw
  """
}
