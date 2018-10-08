#!/usr/bin/env nextflow

// Channel.fromFilePairs("files/import-metadata/rfam/*.{ctl,sql}")
//   .map { it[1] }
//   .set { rfam_files }
// process import_rfam_metadata {
//   echo true
//   input:
//   set file(ctl), file(sql) from rfam_files
//   script:
//   filename = ctl.getName()
//   name = filename.take(filename.lastIndexOf('.'))
//   """
//   set -o pipefail
//   mysql \
//     --host ${params.databases.rfam.mysql.host} \
//     --port ${params.databases.rfam.mysql.port} \
//     --user ${params.databases.rfam.mysql.user} \
//     --database ${params.databases.rfam.mysql.db_name} \
//     < $sql > data.tsv
//     rnac rfam $name data.tsv - | pgloader $ctl
//   """
// }

Channel.fromPath('files/import-metadata/ensembl/metadata/*.sql')
  .map { f -> f.getName().take(f.getName().lastIndexOf('.')) }
  .collectFile(name: "possible-data.txt", newLine: true)
  .set { ensembl_possible }

Channel.empty()
  .mix(Channel.from(params.databases.ensembl.mysql))
  .mix(Channel.from(params.databases.ensembl_genomes.mysql))
  .combine(Channel.fromPath('files/import-metadata/ensembl/analyzed.sql'))
  .combine(ensembl_possible)
  .set { ensembl_info }

process find_ensembl_databases {
  input:
  set val(mysql), file(imported_sql), file(possible) from ensembl_info

  output:
  set val(mysql), file('selected.csv') into ensembl_databases

  """
  psql -f "$imported_sql" "$PGDATABASE" > done.csv

  echo 'show databases' |\
  mysql \
    --host ${mysql.host} \
    --port ${mysql.port} \
    --user ${mysql.user} |\
  rnac ensembl select-tasks - $possible done.csv > selected.csv
  """
}

ensembl_databases
  .flatMap { mysql, csv ->
    data = []
    csv.eachLine { line -> data << ([mysql] + line.tokenize(',')) }
    data
  }
  .map { mysql, name, db -> [mysql, name, db, file("files/import-metadata/ensembl/metadata/${name}.sql")] }
  .set { ensembl_metadata_dbs }

process fetch_internal_ensembl_data {
  maxForks params.import_metadata.ensembl.max_forks

  input:
  set val(mysql), val(name), val(db), file(sql) from ensembl_metadata_dbs

  output:
  set val(name), file('data.csv') into ensembl_metadata
  set val(imported) into ensembl_imports

  script:
  imported = [db, name].join(',')
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
  .map { it -> [file("files/import-metadata/ensembl/metadata/${it[0]}.ctl"), it[1]] }
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
   pgloader --on-error-stop $ctl
   """
}

process mark_imports {
  echo true

  input:
  file('data.txt') from ensembl_imports.collectFile(name: 'data.txt')
  file(ctl) from Channel.fromPath('files/import-metadata/ensembl/mark-analyzed.ctl')

  """
  cp $ctl local_$ctl
  pgloader --on-error-stop local_$ctl
  """
}
