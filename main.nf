#!/usr/bin/env nextflow

def any_database(String... names) { names.any { n -> params.import_data.databases[n] } }

def as_mysql_cmd = { db ->
  rest = db.db_name ? " --database ${db.db_name}" : ''
  "mysql --host ${db.host} --port ${db.port} --user ${db.user} $rest"
}

// ===========================================================================
// Setup for running the pipeline
// ===========================================================================

// This channel will store all the outputs for all data to import. Doing this
// prevents us from having many large duplicated blocks. Also adding a new file
// type to import is much easier this way, as we only have to modify in a few
// places instead of every possible process that could produce it.
raw_output = Channel.empty()

// This contains what remote data sources we have to fetch. It focuses on the
// main data files, not the 'extra' ones. As the extra ones generally need a
// different method for each database and have few commonalities.
to_fetch = Channel.empty()

for (entry in params.import_data.databases) {
  if (!entry.value) {
    continue
  }

  name = entry.key
  database = params.databases[name]

  // The custom property should be a list of directories that contain csv files to
  // import. For example if custom was ./custom we would expect it to contain:
  // ./custom/accessions.csv, ./custom/genomic_locations.csv, etc. The files may be
  // nested so: ./custom/mirbase/accessions.csv would also work.
  if (name == "custom") {
    raw_output.mix(Channel.fromPath("${base}/**.csv")).flatten().set { raw_output }
  } else {
    rnc_cmd = "${database.get('group', 'external')} ${name}"
    db_channel = Channel.value([rnc_cmd, database.fetch])
    to_fetch.mix(db_channel).set { to_fetch }
  }
}

//=============================================================================
// Fetch data as well as all extra data for import
//=============================================================================

process fetch_data {
  input:
  set val(name), val(incomplete) from to_fetch

  output:
  set val(name), file("${database.pattern}") into fetched

  script:
  database = ["excluded_patterns": []] + incomplete
  database.cmd = database.cmd ? database.cmd : "fetch '${database.remote}' '${database.pattern}'"
  clean_cmd = ''
  for (rm_pattern in database['excluded_patterns']) {
    clean_cmd = find_cmd + "find . -name '$rm_pattern}' | xargs rm\n"
  }

  """
  ${database.cmd}
  $clean_cmd
  """
}

process fetch_pdb_extra {
  when:
  any_database('pdb')

  output:
  set val('pdb'), file('pdb-extra.json') into pdb_extra

  """
  rnac pdb extra pdb-extra.json
  """
}

process fetch_rfam_extra {
  when:
  any_database('rfam', 'ensembl')

  input:
  file(sql) from Channel.fromPath("files/import-data/rfam/*.sql")

  output:
  set val(name), file('data.tsv') into all_rfam_extra
  file(data_name) into rfam_output
  set val(name), file(data_name) into rfam_terms

  script:
  filename = ctl.getName()
  name = filename.take(filename.lastIndexOf('.'))
  data_name = "rfam-${name}.csv"
  """
  set -o pipefail

  ${as_mysql_cmd(params.databases.rfam.mysql)} < $sql > data.tsv
  rnac rfam $name data.tsv $data_name
  """
}

all_rfam_extra
  .filter { n, f -> n == 'families' }
  .flatMap { _, f ->
    results = []
      ['rfam', 'ensembl'].each { db ->
        if (params.import_data.databases[db]) {
          results << [db, f]
        }
        results
      }
  }
  .set { rfam_based_extra }

rfam_terms
  .filter { n, f -> n == "ontology-terms" }
  .map { n, f -> f }
  .set { rfam_terms }

Channel
  .from(params.databases.ena.tpa_urls)
  .collectFile(name: "urls", newLine: true)
  .set { tpa_url_file }

process fetch_ena_extra {
  when:
  any_database('ena')

 input:
 file tpa_file from tpa_url_file

 output:
 set val('ena'), file("tpa.tsv") into ena_extra

  """
  cat $tpa_file | xargs wget -O - >> tpa.tsv
  """
}

process fetch_rgd_extra {
  when:
  any_database('rgd')

  output:
  set val('rgd'), file('genes.txt') into rgd_extra

  """
  wget -O genes.txt ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt
  """
}

Channel.fromPath('files/import-data/ensembl/*.sql')
  .map { f -> f.getName().take(f.getName().lastIndexOf('.')) }
  .collectFile(name: "possible-data.txt", newLine: true)
  .set { ensembl_possible }

Channel.from(params.databases.ensembl.mysql, params.databases.ensembl_genomes.mysql)
  .combine(Channel.fromPath('files/import-data/ensembl-analyzed.sql'))
  .combine(ensembl_possible)
  .set { ensembl_info }

process find_ensembl_tasks {
  when:
  any_database('ensembl')

  input:
  set val(mysql), file(imported_sql), file(possible) from ensembl_info

  output:
  set val(mysql), file('selected.csv') into ensembl_databases

  """
  set -o pipefail

  psql -v ON_ERROR_STOP=1 -f "$imported_sql" "$PGDATABASE" > done.csv
  echo 'show databases' | ${as_mysql_cmd(mysql)} > dbs.txt
  rnac ensembl select-tasks dbs.txt $possible done.csv > selected.csv
  """
}

Channel.empty()
  .mix(ena_extra)
  .mix(rfam_based_extra)
  .mix(rgd_extra)
  .mix(pdb_extra)
  .groupBy()
  .set { extra }

fetched
  .flatMap { name, filenames ->
    filenames = filenames.class.isArray() ? filenames : [filenames]
    results = []
    filenames.each { filename ->
      results << [name, filename]
    }
    results
  }
  .combine(extra)
  .map { name, data_file, extra ->
    // Pretty sure this is because groovy is messing with types or something
    to_add = extra[name] ? extra[name][0][1..-1] : []
    [name, data_file, to_add]
  }
  .set { to_process }

ensembl_databases
  .flatMap { mysql, csv ->
    data = []
    csv.eachLine { line -> data << ([mysql] + line.tokenize(',')) }
    data
  }
  .map { mysql, name, db -> [mysql, name, db, file("files/import-data/ensembl/${name}.sql")] }
  .set { ensembl_metadata_dbs }

//=============================================================================
// Process data
//=============================================================================

process process_data {
  memory params.databases[name].get('memory', '2 GB')

  input:
  set val(name), file(input_file), file(extra) from to_process

  output:
  file "*.csv" into all_processed_output mode flatten

  script:
  if (input_file.toString().endsWith('.gz')) {
    """
    set -o pipefail

    zcat $input_file | rnac ${name} - ${extra}
    """
  } else {
    """
    rnac ${name} ${input_file} ${extra}
    """
  }
}

all_processed_output
  .view()
  .into { processed_output; all_terms }

all_terms
  .filter { f -> f.getName() == "terms.csv" }
  .mix(rfam_terms)
  .collect()
  .set { terms }

process fetch_ontology_information {
  input:
  file('terms*.csv') from terms

  output:
  file('ontology_terms.csv') into term_info

  script:
  """
  sort -u terms*.csv > unique-terms.txt
  rnac ontologies lookup-terms unique-terms.txt ontology_terms.csv
  """
}

process process_ensembl_metadata {
  maxForks params.import_metadata.ensembl.max_forks

  input:
  set val(mysql), val(name), val(db), file(sql) from ensembl_metadata_dbs

  output:
  file("${name}.csv") optional true into ensembl_output
  val("$name,$db") into ensembl_imported

  script:
  """
  ${as_mysql_cmd(mysql)} -N --database $db < $sql > data.tsv
  rnac ensembl $name data.tsv ${name}.csv
  """
}

ensembl_imported
  .collectFile(name: "ensembl_imported.csv", newLine: true)
  .set { ensembl_imported_output }

//=============================================================================
// Import and release data
//=============================================================================

raw_output
  .mix(processed_output)
  .mix(rfam_output)
  .mix(term_info)
  .mix(ensembl_output)
  .mix(ensembl_imported_output)
  .map { f ->
    filename = f.getName()
    name = filename.take(filename.lastIndexOf('.'))
    ctl = file("files/import-data/load/${name.replace('_', '-')}.ctl")
    [[name, ctl], f]
  }
  .filter { it[0][1].exists() }
  .groupTuple()
  .map { it -> [it[0][0], it[0][1], it[1]] }
  .set { to_load }

to_load.println()

// process merge_and_import {
//   echo true

//   input:
//   set val(name), file(ctl), file('raw*.csv') from to_load

//   output:
//   val(name) into loaded

//   """
//   split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
//   """
// }

// loaded
//   .flatMap { name -> file("files/import-data/post-release/*__${name.replace('_', '-')}.sql") }
//   .mix(Channel.fromPath('files/import-data/post-release/99__cleanup.sql'))
//   .filter { f -> f.exists() }
//   .collect()
//   .set { post_scripts }

// process release {
//   echo true
//   maxForks 1

//   input:
//   file(post) from post_scripts

//   output:
//   val('done') into post_release

//   """
//   set -o pipefail

//   rnac run-release
//   find . -name '*.sql' -print0 |\
//   sort -z |\
//   xargs -r0 -I {} psql -v ON_ERROR_STOP=1 -f {} "$PGDATABASE"
//   """
// }

// post_release
//   .ifEmpty('no release')
//   .set { sequences_to_scan }
