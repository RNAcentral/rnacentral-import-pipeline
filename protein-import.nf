#!/usr/bin/env nextflow

process find_ensembl_databases {
  output:
  stdout into ensembl_databases

  """
  echo 'show databases' |\
  mysql \
    --host ${params.databases.ensembl.mysql.host} \
    --port ${params.databases.ensembl.mysql.port} \
    --user ${params.databases.ensembl.mysql.user} |\
  max-ensembl-database
  """
}

ensembl_databases
  .splitCsv()
  .combine(Channel.fromPath('files/protein-info/ensembl.sql'))
  .set { ensembl_data }

process ensembl_protein_info {
  maxForks params.protein_import.ensembl.max_forks

  input:
  set val(db), file(sql) from ensembl_data

  output:
  file('proteins.tsv') into ensembl_proteins

  """
  mysql -N \
    --host ${params.databases.ensembl.mysql.host} \
    --port ${params.databases.ensembl.mysql.port} \
    --user ${params.databases.ensembl.mysql.user} \
    --database ${db} \
    < $sql > proteins.tsv
  """
}

process import_ensembl_proteins {
   input:
   file('proteins*.tsv') from ensembl_proteins.collect()
   file(ctl) from Channel.fromPath('files/protein-info/ensembl.ctl')

   """
   cat proteins*.tsv | rnac proteins ensembl - - | pgloader $ctl
   """
}
