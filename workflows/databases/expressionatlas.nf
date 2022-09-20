include { slack_message } from './workflows/utils/slack'
include { slack_closure } from './workflows/utils/slack'

process fetch_data {
  queue datamover

  output:
  path('*.tsv')

  script:
  """
  BASE_DIR="/nfs/ftp/public/databases/microarray/data/atlas/experiments/"
  find $BASE_DIR -type f -iname "*tpms.tsv" -or -iname "*analytics.tsv" -or -iname "*sdrf.tsv" > all_relevant_files

  while read f; do
    cp $f .
    echo "$f"
  done < all_relevant_files
  """
}

process fetch_lookup {
  queue short

  input:
    path(query)

  output:
    path("lookup_dump.csv")

    """
    psql -f $query $PGDATABASE > lookup_dump.csv
    """
}


process parse_tsvs {
  memory 50.GB

  input:
  path(tsvs), path(lookup)

  output:
  path('*.csv')

  """
  expression-parse -i $tsvs -o all_genes.csv -s GeneID experiment
  expression-parse lookup -g all_genes.csv -l $lookup -o exp_parse_stage2.json
  rnac expressionatlas parse exp_parse_stage2.json .
  """

}


workflow expressionatlas {

  emit: data
  main:

    Channel.of("Expression Atlas import starting...") | slack_message
    Channel.fromPath('files/import-data/expressionatlas/lookup-dump-query.sql') | set { lookup_sql }

    lookup_sql | fetch_lookup | set { lookup }
    fetch_data | set { tsvs }

    tsvs| combine(lookup) | parse_tsvs | set { data }

}
