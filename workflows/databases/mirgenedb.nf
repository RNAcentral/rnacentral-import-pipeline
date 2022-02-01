process mirgenedb {
  when { params.databases.mirgenedb.run }

  output:
  path('*.csv')

  """
  scp $params.databases.mirgenedb.remote mirgenedb.json
  psql \
    --csv \
    --command='select assembly_id,assembly_ucsc from ensembl_assembly where assembly_ucsc is not null'
    "$PGDATABASE" > assemblies.tsv
  rnac mirgenedb parse assemblies.tsv mirgenedb.json .
  """
}
