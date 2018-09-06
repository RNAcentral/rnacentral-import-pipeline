#!/usr/bin/env nextflow

Channel.fromFilePairs("files/import-metadata/rfam/*.{ctl,sql}")
  .map { it[1] }
  .set { rfam_files }

process import_rfam_metadata {
  input:
  set file(ctl), file(sql) from rfam_files

  """
  mysql \
    --host ${params.databases.rfam.mysql.host} \
    --port ${params.databases.rfam.mysql.port} \
    --user ${params.databases.rfam.mysql.user} \
    --database ${params.databases.rfam.mysql.db_name} \
    < $sql | pgloader $ctl
  """
}

process fetch_crs {
  input:
  val(remote) from Channel.from(params.crs.path)

  output:
  file('*.tsv.gz') into raw_crs mode flatten

  script:
  """
  fetch "$remote" "*.tsv.gz"
  """
}

process process_crs {
  input:
  file(crs) from raw_crs

  output:
  file('complete_features.csv') into processed_crs

  """
  zcat $crs | rnac extra crs - complete_features.csv
  """
}

process import_crs {
  input:
  file('complete_features*.csv') from processed_crs.collect()
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
