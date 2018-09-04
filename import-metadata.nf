#!/usr/bin/env nextflow

process import_rfam_metadata {
  input:
  set file(ctl), file(sql) from Channel.fromFilePairs("files/import-data/rfam/*.{ctl,sql}")

  """
  mysql \
    --host ${params.databases.rfam.mysql.host} \
    --port ${params.databases.rfam.mysql.port} \
    --user ${params.databases.rfam.mysql.user} \
    --database ${params.databases.rfam.mysql.db_name} \
    < $sql | pgloader $ctl
  """
}

process import_crs {
  input:
  file(crs) from Channel.fromPath(params.crs.path)

  output:
  file('complete_features.csv') into raw_crs

  """
  zcat $crs | rnc extra crs - complete_features.csv
  """
}

process import_crs {
  input:
  file('complete_features*.csv') from raw_crs.collect()
  file(ctl) from Channel.fromPath('files/import-metadata/crs-features.ctl')

  """
  cp $ctl crs.ctl
  pgloader crs.ctl
  """
}

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
   sort proteins*.tsv | rnac proteins ensembl - - | pgloader $ctl
   """
}
