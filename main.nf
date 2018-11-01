#!/usr/bin/env nextflow

def any_database(String... names) { names.any { n -> params.import_data.databases[n] } }

def as_mysql_cmd = { db ->
  rest = db.db_name ? " --database ${db.db_name}" : ''
  "mysql --host ${db.host} --port ${db.port} --user ${db.user} $rest"
}

// ===========================================================================
// Compute initial tasks
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

// If we are processing Rfam, Ensembl, or GENCODE we must fetch and import some
// Rfam based metadata. These databases use Rfam in one way or another so we
// should stay up-to-date with the current Rfam data. For example they may be
// creating annotations of a partial lncRNA which we do not want to import. By
// tracking the latest Rfam we will know what is a partial lncRNA. This is done
// by querying the Rfam MySQL database and then using the Rfam commands to
// process the query results. We can then import the processed data as normal.
if (any_database('rfam', 'ensembl', 'gencode')) {
  Channel.fromPath("files/import-data/rfam/*.sql")
    .map { query ->
      name = query.getBaseName()
      fetch = rfam_mysql + [cmd: 'mysql', query: query, pattern: 'data.tsv']
      ["rfam $name", fetch]
    }
    .set { rfam_queries }
  to_fetch.mix(rfam_queries).set { to_fetch }
}

Channel.fromPath('files/import-data/ensembl/*.sql')
  .map { f -> f.getBaseName() }
  .collectFile(name: "possible-data.txt", newLine: true)
  .set { ensembl_possible }

Channel.from(params.databases.ensembl.mysql, params.databases.ensembl_genomes.mysql)
  .combine(Channel.fromPath('files/import-data/ensembl-analyzed.sql'))
  .combine(ensembl_possible)
  .set { ensembl_info }

process find_ensembl_tasks {
  input:
  set val(mysql), file(imported_sql), file(possible) from ensembl_info

  output:
  set val(mysql), file('selected.csv') into ensembl_task_summary

  """
  set -o pipefail

  psql -v ON_ERROR_STOP=1 -f "$imported_sql" "$PGDATABASE" > done.csv
  echo 'show databases' | ${as_mysql_cmd(mysql)} > dbs.txt
  rnac ensembl select-tasks dbs.txt $possible done.csv > selected.csv
  """
}

//=============================================================================
// Fetch data as well as all extra data for import
//=============================================================================

process fetch_data {
  input:
  set val(name), val(incomplete) from to_fetch

  output:
  set val(name), file("${database.pattern}") into all_fetched

  script:
  database = ["excluded_patterns": []] + incomplete
  if (database.cmd == "mysql") {
    """
    ${as_mysql_cmd(database)} < ${database.query} > ${database.pattern}
    """
  } else {
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
}

all_fetched
  .into { fetched; for_gencode; rfam_based }

process fetch_pdb_extra {
  when:
  any_database('pdb')

  output:
  set val('external pdb'), file('pdb-extra.json') into pdb_extra

  """
  rnac pdb extra pdb-extra.json
  """
}

// This filters the fetched data to extract the Rfam families data. This data is
// used by the Ensembl and GENCODE processing logic later on. Thus we create a
// channel of just that data after fetching it.
rfam_based
  .filter { n, f -> n == "rfam families" }
  .flatMap { n, f ->
    results = []
      ['ensembl', 'gencode'].each { db ->
        if (params.import_data.databases[db]) {
          results << [db, f]
        }
        results
      }
  }
  .set { rfam_based_extra }

// Here we split the fetched data from ensembl to pull out all Ensembl human and
// mouse raw data. This data is parsed both for Ensembl import and for GENCODE
// import. The raw GENCODE GFF3 file doesn't contain all useful information like
// sequences and names.
for_gencode
  .filter { n, fs -> any_database('gencode') && n == "external ensembl" }
  .flatMap { n, fs -> fs }
  .filter { f ->
    species = f.getBaseName()
    species.startsWith('Homo_sapiens') || species.startsWith('Mus_mus')
  }
  .map { f -> ['external gencode', f] }
  .set { gencode_extra }

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
 set val('external ena'), file("tpa.tsv") into ena_extra

  """
  cat $tpa_file | xargs wget -O - >> tpa.tsv
  """
}

process fetch_rgd_extra {
  when:
  any_database('rgd')

  output:
  set val('external rgd'), file('genes.txt') into rgd_extra

  """
  wget -O genes.txt ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt
  """
}

//=============================================================================
// Create channels for all data processing tasks
//=============================================================================

ensembl_task_summary
  .flatMap { mysql, csv ->
    data = []
    csv.eachLine { line ->
      (name, db) = line.tokenize(',')
      updated = mysql + [db_name: db]
      data << [updated, name, file("files/import-data/ensembl/${name}.sql")]
    }
    data
  }
  .set { ensembl_metadata_tasks }

Channel.empty()
  .mix(ena_extra)
  .mix(gencode_extra)
  .mix(rfam_based_extra)
  .mix(rgd_extra)
  .mix(pdb_extra)
  .groupBy()
  .set { extra }

fetched
  .combine(extra)
  .flatMap { name, filenames, extra ->
    // Pretty sure this weirdness is because groovy is messing with types or something
    to_add = extra[name] ? extra[name][0][1..-1] : []
    fs = [filenames].flatten()
    results = []
    fs.each { f -> results << [name, f, to_add] }
    results
  }
  .set { to_process }

//=============================================================================
// Process data
//=============================================================================

process process_data {
  memory params.databases[name].get('memory', '2 GB')

  input:
  set val(name), file(input_file), file(extra) from to_process

  output:
  file "*.csv" optional true into all_processed_output mode flatten

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

processed_output = Channel.create()
terms = Channel.create()

all_processed_output
  .choice(terms, refs, processed_output) { f ->
    names = ["terms.csv", "refs.csv"]
    index = names.indexOf(f.getName())
    return index >= 0 : index : names.size()
  }

refs.collect()
  .set { refs_to_process }

process fetch_ontology_information {
  input:
  file('terms*.csv') from terms.collect()

  output:
  file('ontology_terms.csv') into term_info

  script:
  """
  sort -u terms*.csv > unique-terms.txt
  rnac ontologies lookup-terms unique-terms.txt ontology_terms.csv
  """
}

process fetch_publications {
  intput:
  set file('ref_ids*.csv') from refs_to_process

  output:
  file('references.csv') into references_output

  script:
  pubs = params.import_metadata.publications
  """
  set -o pipeline

  find . -name 'ref_ids*.csv' | xargs cat >> all-ids
  cut -d, -f1 all-ids | uniq | sort -u > ref_ids

  lines="$(wc -l ref_ids | cut -d ' ' -f1)"
  rnac publications fetch all-ids references.csv
  """
}

process process_ensembl_metadata {
  maxForks params.import_metadata.ensembl.max_forks

  input:
  set val(mysql), val(name), file(sql) from ensembl_metadata_tasks

  output:
  file("${name}.csv") optional true into ensembl_output
  val("$name,$db") into ensembl_imported

  script:
  """
  ${as_mysql_cmd(mysql)} -N < $sql > data.tsv
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
  // .mix(references_output)
  .mix(ensembl_output)
  .mix(ensembl_imported_output)
  .map { f ->
    name = f.getBaseName()
    ctl = file("files/import-data/load/${name.replace('_', '-')}.ctl")
    [[name, ctl], f]
  }
  .filter { it[0][1].exists() }
  .groupTuple()
  .map { it -> [it[0][0], it[0][1], it[1]] }
  .set { to_load }

process merge_and_import {
  echo true

  input:
  set val(name), file(ctl), file('raw*.csv') from to_load

  output:
  val(name) into loaded

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  """
}

loaded
  .flatMap { name -> file("files/import-data/post-release/*__${name.replace('_', '-')}.sql") }
  .mix(Channel.fromPath('files/import-data/post-release/99__cleanup.sql'))
  .filter { f -> f.exists() }
  .collect()
  .set { post_scripts }

process release {
  echo true
  maxForks 1

  input:
  file(post) from post_scripts

  output:
  val('done') into post_release

  """
  set -o pipefail

  rnac run-release
  find . -name '*.sql' -print0 |\
  sort -z |\
  xargs -r0 -I {} psql -v ON_ERROR_STOP=1 -f {} "$PGDATABASE"
  """
}

post_release
  .ifEmpty('no release')
  .into { flag_for_qa; flag_for_mapping }

//=============================================================================
// QA scans
//=============================================================================

flag_for_qa
  .combine(Channel.fromPath('files/qa/rfam-scan.sql'))
  .set { qa_queries }

process fetch_sequences {
  input:
  set val(status), file(query) from qa_queries

  output:
  file('parts/*.fasta') into sequences_to_scan mode flatten

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${params.qa.rfam_scan.chunk_size} --out-dir 'parts/' shuffled.fasta
  """
}

sequences_to_scan
  .combine(Channel.fromPath(params.qa.rfam_scan.cm_files).collect())
  .set { sequences_to_scan }

process infernal_scan {
  queue 'mpi-rh7'
  cpus params.qa.rfam_scan.cpus
  clusterOptions "-M ${params.qa.rfam_scan.cm_memory} -R 'rusage[mem=${params.qa.rfam_scan.cm_memory}]' -a openmpi"
  module 'mpi/openmpi-x86_64'

  input:
  set file('sequences.fasta'), file(cm_file) from sequences_to_scan

  output:
  set val('rfam'), file('hits.csv') into processed_rfam_hits

  """
  mpiexec -mca btl ^openbib -np ${params.qa.rfam_scan.cpus} \
  cmscan \
    -o output.inf \
    --tblout results.tblout \
    --clanin ${params.qa.rfam_scan.clans} \
    --oclan \
    --fmt 2 \
    --acc \
    --cut_ga \
    --rfam \
    --notextw \
    --nohmmonly \
    --mpi \
    "Rfam.cm" \
    sequences.fasta
  rnac qa tblout2csv results.tblout hits.csv
  """
}

processed_rfam_hits
  .groupTuple()
  .combine(Channel.fromPath('files/qa/rfam-scan.ctl'))
  .set { hits_to_import }

process import_qa_data {
  echo true

  input:
  set val(name), file('raw*.csv'), file(ctl) from hits_to_import

  output:
  val('done') into qa_imported

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  """
}

qa_imported
  .ifEmpty('no qa')
  .set { qa_imported }
