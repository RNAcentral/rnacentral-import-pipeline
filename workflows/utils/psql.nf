process simple_query {
  input:
  file(query)

  output:
  file('result')

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > result
  """
}
