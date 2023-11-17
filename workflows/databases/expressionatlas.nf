process fetch_data {
  queue 'datamover'
  container ''
  errorStrategy 'ignore'
  cpus 10

  input:
    val base_dir

  output:
  path('tsv_files')

  """
  mkdir tsv_files
  find "${base_dir}" -type f -name "*.tsv" -print0 | xargs -0 -I {} -P 20 cp {} tsv_files || true
  find "${base_dir}" -type f -name "*configuration.xml" -print0 | xargs -0 -I {} -P 20 cp {} tsv_files || true
  """
}

process fetch_lookup {

  input:
    path (query)

  output:
    path("lookup_dump.csv")

    """
    psql -f $query $PGDATABASE > lookup_dump.csv
    """
}


process parse_tsvs {
  memory { 64.GB * task.attempt }
  cpus 16
  container ''
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

  input:
  path(tsvs)

  output:
  path('chunk_*')

  """
  expression-parse parse -i $tsvs -o all_genes.csv
  parallel --block 500M -a all_genes.csv --header : --pipepart 'cat > chunk_{#}'
  """

}

process lookup_genes {
  memory { 32.GB * task.attempt }
  cpus 16
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  input:
    tuple path(lookup), path(genes)

  output:
    path('*.csv')

  """
  expression-parse lookup -g $genes -l $lookup -o exp_parse_stage2.json
  rnac expressionatlas parse exp_parse_stage2.json .
  """
}


workflow expressionatlas {

  emit: data
  main:

  if( params.databases.expressionatlas.run ) {
    Channel.fromPath('files/import-data/expressionatlas/lookup-dump-query.sql') | set { lookup_sql }
    Channel.of(params.databases.expressionatlas.remote) | set { tsv_path }
    lookup_sql | fetch_lookup | set { lookup }
    tsv_path \
    | fetch_data \
    | parse_tsvs \
    | flatten | set { genes }

   lookup.combine(genes) | lookup_genes \
   | collect \
   | set { data }

  }
  else {
    Channel.empty() | set { data }
  }

}
