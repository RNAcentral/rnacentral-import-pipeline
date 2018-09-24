#!/usr/bin/env nextflow

Channel.fromFilePairs("files/import-metadata/rfam/*.{ctl,sql}")
  .map { it[1] }
  .set { rfam_files }

process import_rfam_metadata {
  echo true

  input:
  set file(ctl), file(sql) from rfam_files

  script:
  filename = ctl.getName()
  name = filename.take(filename.lastIndexOf('.'))
  """
  set -o pipeline

  mysql \
    --host ${params.databases.rfam.mysql.host} \
    --port ${params.databases.rfam.mysql.port} \
    --user ${params.databases.rfam.mysql.user} \
    --database ${params.databases.rfam.mysql.db_name} \
    < $sql > data.tsv
    rnac rfam $name data.tsv - | pglaoder $ctl
  """
}

Channel.empty()
  .mix(Channel.from(params.databases.ensembl.mysql))
  .mix(Channel.from(params.databases.ensembl_genomes.mysql))
  .set { ensembl_info }

process find_ensembl_databases {
  input:
  val(mysql) from ensembl_info

  output:
  set val(mysql), file('selected.csv') into ensembl_databases

  """
  echo 'show databases' |\
  mysql \
    --host ${mysql.host} \
    --port ${mysql.port} \
    --user ${mysql.user} |\
  max-ensembl-database > selected.csv
  """
}

ensembl_databases
  .flatMap { mysql, csv ->
    data = []
    csv.eachLine { line -> data << ([mysql] + line.tokenize(',')) }
    data
  }
  .combine(Channel.fromPath('files/import-metadata/ensembl/*.sql'))
  .map { mysql, db, sql ->
    filename = sql.getName()
    name = filename.take(filename.lastIndexOf('.'))
    [mysql, name, db, sql]
  }
  .set { ensembl_metadata_dbs }

process ensembl_protein_info {
  maxForks params.import_metadata.ensembl.max_forks

  input:
  set val(mysql), val(name), val(db), file(sql) from ensembl_metadata_dbs

  output:
  set val(name), file('data.csv') into ensembl_metadata

  """
  mysql -N \
    --host ${mysql.host} \
    --port ${mysql.port} \
    --user ${mysql.user} \
    --database ${db} \
    < $sql > data.tsv
  rnac ensembl $name data.tsv data.csv
  """
}

ensembl_metadata
  .groupTuple()
  .map { it -> [file("files/import-metadata/ensembl/${it[0]}.ctl"), it[1]] }
  .set { ensembl_loadable }

process import_ensembl_data {
   echo true

   input:
   set file(ctl), file('data*.csv') from ensembl_loadable

   """
   sort -u data*.csv > merged.csv
   cp $ctl local_$ctl
   pgloader local_$ctl
   """
}
