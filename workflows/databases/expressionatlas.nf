include { slack_message } from './workflows/utils/slack'
include { slack_closure } from './workflows/utils/slack'

process fetch_data {
  queue datamover

  output:
  path('*.tsv')

  script:
  """
  BASE_DIR="/nfs/ftp/public/databases/microarray/data/atlas/experiments/"
  find $BASE_DIR -type f -iname "*tpms.tsv" -or -iname "*precentile-ranks.tsv" -or -iname "*sdrf.tsv" > all_relevant_files

  while read f; do
    cp $f .
    echo "$f"
  done < all_relevant_files
  """
}


process parse_tsvs {
  memory 50.GB

  input:
  path(tsvs)

  output:
  path('*.csv')

  """
  expression-parse -i $tsvs -o all_genes.csv -s GeneID experiment

  rnac expressionatlas parse all_genes.csv .
  """

}


workflow expressionatlas {

  emit: data
  main:

    Channel.of("Expression Atlas import starting...") | slack_message

    fetch_data | parse_tsvs | set { data }

}
