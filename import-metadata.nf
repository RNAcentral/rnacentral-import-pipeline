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
  set -o pipefail

  mysql \
    --host ${params.databases.rfam.mysql.host} \
    --port ${params.databases.rfam.mysql.port} \
    --user ${params.databases.rfam.mysql.user} \
    --database ${params.databases.rfam.mysql.db_name} \
    < $sql > data.tsv
    rnac rfam $name data.tsv - | pgloader $ctl
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
  rnc ensembl select-databases - > selected.csv
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

process fetch_internal_ensembl_data {
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
   mkdir merged
   find . -name 'data*.csv' |\
   xargs cat |\
   split --additional-suffix=.csv -dC ${params.import_data.chunk_size} - merged/data

   cp $ctl merged/$ctl
   cd merged
   pgloader $ctl
   """
}
