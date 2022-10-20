process fetch_data {
  queue 'datamover'
  container ''
  errorStrategy 'ignore'

  input:
    path("base_dir")

  output:
  path('tsv_files')

  """
  mkdir tsv_files
  find $base_dir -type f .. | grep -v E-MTAB-4748 | xargs -I {} -P 10 cp {} tsv_files
  """
}

process fetch_lookup {
  queue 'short'

  input:
    path (query)

  output:
    path("lookup_dump.csv")

    """
    psql -f $query $PGDATABASE > lookup_dump.csv
    """
}


process parse_tsvs {
  memory 92.GB

  input:
  path(tsvs)
  path(lookup)

  output:
  path('*.csv')

  """
  expression-parse parse -i $tsvs -o all_genes.csv
  wc -l all_genes.csv
  expression-parse lookup -g all_genes.csv -l $lookup -o exp_parse_stage2.json
  rnac expressionatlas parse exp_parse_stage2.json .
  """

}


workflow expressionatlas {

  emit: data
  main:

  if( params.databases.expressionatlas.run ) {
    Channel.fromPath('files/import-data/expressionatlas/lookup-dump-query.sql') | set { lookup_sql }
    Channel.fromPath($params.databases.expressionatlas.remote) | set { tsv_path }
    lookup_sql | fetch_lookup | set { lookup }
    tsv_path | fetch_data | set { tsvs }

   parse_tsvs(tsvs, lookup) | set { data }
  }
  else {
    Channel.empty() | set { data }
  }

}
