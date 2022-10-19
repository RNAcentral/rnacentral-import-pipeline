process fetch_data {
  queue 'datamover'
  container ''
  errorStrategy 'ignore'

  output:
  path('tsv_files')

  """
  find /nfs/ftp/public/databases/microarray/data/atlas/experiments -type f -name '*tpms.tsv' -or -name '*analytics.tsv' -or -name '*sdrf.tsv' -or -name '*configuration.xml' > all_relevant_files || true

  sed -i '/E-MTAB-4748/d' all_relevant_files
  sed -i '/E-GEOD-30281/d' all_relevant_files


  mkdir tsv_files
  cat all_relevant_files | while read line
  do
    cp \$line ./tsv_files
    echo \$line
  done

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

    lookup_sql | fetch_lookup | set { lookup }
    fetch_data | set { tsvs }

   parse_tsvs(tsvs, lookup) | set { data }
  }
  else {
    Channel.empty() | set { data }
  }

}
