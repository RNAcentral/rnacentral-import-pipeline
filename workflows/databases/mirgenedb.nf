process mirgenedb {
  memory { 2.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  when { params.databases.mirgenedb.run }

  output:
  path('*.csv')

  """
  scp $params.databases.mirgenedb.remote mirgenedb.json
  psql \
    --command='COPY (select assembly_id,assembly_ucsc from ensembl_assembly where assembly_ucsc is not null) TO STDOUT (FORMAT CSV)' \
    "$PGDATABASE" > assemblies.tsv
  rnac mirgenedb parse assemblies.tsv mirgenedb.json .
  """
}
