#!/usr/bin/env nextflow

Channel.fromFilePairs("files/import-metadata/rfam/*.{ctl,sql}")
  .map { it[1] }
  .set { rfam_files }

// process import_rfam_metadata {
//   echo true

//   input:
//   set file(ctl), file(sql) from rfam_files

//   """
//   mysql \
//     --host ${params.databases.rfam.mysql.host} \
//     --port ${params.databases.rfam.mysql.port} \
//     --user ${params.databases.rfam.mysql.user} \
//     --database ${params.databases.rfam.mysql.db_name} \
//     < $sql | pgloader $ctl
//   """
// }

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
  .combine(Channel.fromPath('files/import-metadata/ensembl/*.sql'))
  .map { it ->
    filename = it[1].getName()
    name = filename.take(filename.lastIndexOf('.'))
    [name, it[0], it[1]]
  }
  .set { ensembl_metadata_dbs }

process ensembl_protein_info {
  maxForks params.import_metadata.ensembl.max_forks

  input:
  set val(name), val(db), file(sql) from ensembl_metadata_dbs

  output:
  set val(name), file('data.csv') into ensembl_metadata

  """
  mysql -N \
    --host ${params.databases.ensembl.mysql.host} \
    --port ${params.databases.ensembl.mysql.port} \
    --user ${params.databases.ensembl.mysql.user} \
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
