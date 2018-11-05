#!/usr/bin/env nextflow

def any_database(String... names) { names.any { n -> params.import_data.databases[n] } }

def as_mysql_cmd = { db ->
  rest = db.db_name ? " --database ${db.db_name}" : ''
  "mysql --host ${db.host} --port ${db.port} --user ${db.user} $rest"
}

def fetch_for(incomplete) {
  if (!incomplete) {
    return ''
  }

  db = ["cmd": "fetch"] + incomplete
  if (db.cmd == "mysql") {
    cmd = """
    ${as_mysql_cmd(db)} < ${db.query} > ${db.pattern}
    """
  } else {
    ps = database.get('excluded_patterns', [])
    clean = ps.inject('', { agg, p -> agg + "find . -name '$p' | xargs rm\n" })
    cmd = """\
    ${db.cmd} ${db.name.join(' ')} '${db.remote}' '${db.pattern}'
    ${clean}
    """
  }
  return cmd.stripIndent()
}

def process_using = { db ->
  input_file = file(db.input_file)
  name = db.name.join(' ')
  extra = db.extra ? db.extra : ''
  if (input_file.getName().endsWith('.gz')) {
    return "zcat ${db.input_file} | rnac ${name} - ${extra}"
  }
  return "rnac ${name} ${db.input_file} ${extra}"
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
// main data files, not the 'extra' ones. Separating the two makes rejoining
// them later much easier.
to_fetch = Channel.empty()

// Here we put the 'extra' data files to fetch.
to_fetch_extra = Channel.empty()

// This contains the tasks for jobs that can be fetched and processed in a
// single task. Some databases, like all JSON based import ones provide a single
// file that we need to process. In that case we don't have to split the
// fetching and processing into separate tasks. We can do 1 task instead of 2
// (makes the job scheduler and ITS happier) and perform a single task that does
// both.
to_fetch_and_process = Channel.empty()

// GENCODE relies upon Ensembl so ensure that Ensembl is always processed if
// GENCODE is.
if (params.import_data.databases.gencode) {
  params.import_data.datatabases.ensembl = true
}

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
    continue
  }

  fetch = [name: [database.get('group', 'external'), name]] + database.fetch

  extra = null
  if (database.extra) {
    extra = [name: ['extra', name]] + database.extra
    if (extra.remotes) {
      extra.remote = file("$name-remotes.txt")
      extra.remote.text = extra.remotes.join('\n')
    }
  }

  // Here we detect if a database only requires parsing a single file. This is
  // done by ensure there are no extra files and there is no '*' in the output
  // pattern name. So long as that is true we can then place this in a channel
  // for fetching and processing as a single task. This also checks if the
  // as_single_task flag is set to false. This is something that GENCODE has
  // to set because it also uses Ensembl data. However, we don't add any
  // extras to the configuration because then we would end up fetching things
  // twice. We can just reuse what we fetch as seen below.
  no_merge = database.fetch.pattern.contains('*') || !database.get('as_single_task', true)
  if (no_merge) {
    to_fetch.mix(Channel.value(fetch)).set { to_fetch }
    to_fetch_extra = extra ? to_fetch_extra.mix(Channel.value(extra)) : to_fetch_extra

  // If we can merge the fetch and process steps into a single step then we do
  // so here.
  } else {
    fetch.extra = extra
    to_fetch_and_process
      .mix(Channel.value(fetch))
      .set { to_fetch_and_process }
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
  executor 'local'

  when:
  any_database('ensembl')

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
  tag name

  input:
  val(db) from to_fetch

  output:
  set val(name), file("${db.pattern}") into all_fetched

  script:
  name = db.name.join(' ')
  """
  ${fetch_for(db)}
  """
}

all_fetched.into { fetched; rfam_based; for_gencode }

process fetch_extra_data {
  tag name

  input:
  val(db) from to_fetch_extra

  output:
  set val(name), file("${db.pattern}") into fetched_extra

  script:
  name = "external ${db.name[1]}"
  """
  ${fetch_for(db)}
  """
}

// This filters the fetched data to extract the Rfam families data. This data is
// used by the Ensembl and GENCODE processing logic later on. Thus we create a
// channel for that data after fetching it, and make it into an extra set of
// data.
rfam_based
  .filter { n, f -> n == "rfam families" }
  .flatMap { n, f ->
    dbs = params.import_data.databases
    ['ensembl', 'gencode'].inject([], { res, db -> dbs[db] ? res + [db, f] : res })
  }
  .set { rfam_based_extra }

// Here we split the fetched data from Ensembl to pull out all Ensembl human and
// mouse raw data. This data is reused for GENCODE import. We need the EMBL
// files as an extra file to extract all useful information like names and such.
for_gencode
  .filter { n, fs -> any_database('gencode') && n == "external ensembl" }
  .flatMap { n, fs -> fs }
  .filter { f ->
    species = f.getBaseName()
    species.startsWith('Homo_sapiens') || species.startsWith('Mus_mus')
  }
  .map { f -> ['external gencode', f] }
  .set { gencode_extra }

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
  .mix(fetched_extra)
  .mix(gencode_extra)
  .mix(rfam_based_extra)
  .groupBy()
  .set { extra }

fetched
  .combine(extra)
  .flatMap { name, filenames, extra ->
    // Pretty sure this weirdness is because groovy is messing with types or something
    to_add = extra[name] ? extra[name][0][1..-1] : []
     [filenames].flatten().inject([], { agg, f -> agg << [name, f, to_add] })
  }
  .view { "Will process: ${it}" }
  .set { to_process }

//=============================================================================
// Process data
//=============================================================================

process fetch_and_process {
  tag name
  memory { params.databases[incomplete.name[1]].get('memory', '2 GB') }

  input:
  val(incomplete) from to_fetch_and_process

  output:
  file "*.csv" into all_fetched_and_processed_output mode flatten

  script:
  db = [input_file: incomplete.pattern] + incomplete
  """
  set -o pipefail

  ${fetch_for(db)}
  ${fetch_for(db.extra)}
  ${process_using(db)}
  """
}

process process_data {
  tag name
  memory { params.databases[incomplete.name[1]].get('memory', '2 GB') }

  input:
  set val(name), file(input_file), file(extra) from to_process

  output:
  file "*.csv" into all_processed_output mode flatten

  script:
  db = [name: name, input_file: input_file, extra: extra]
  """
  set -o pipefail

  ${process_using(db)}
  """
}

processed_output = Channel.create()
refs = Channel.create()
terms = Channel.create()

all_processed_output
  .mix(all_fetched_and_processed_output)
  .choice(terms, refs, processed_output) { f ->
    names = ["terms.csv", "ref_ids.csv"]
    index = names.indexOf(f.getName())
    return index >= 0 ? index : names.size()
  }

process lookup_ontology_information {
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

process lookup_publications {
  input:
  file('ref_ids*.csv') from refs.collect()

  output:
  file('references.csv') into references_output

  script:
  """
  set -o pipeline

  find . -name 'ref_ids*.csv' | xargs cat | sort -u >> all-ids

  rnac publications fetch all-ids references.csv
  """
}

// This must be it's own task becuase we need to have a specific number of forks
// and this cannot be set dynamically. Otherwise we could just put it into the
// generic fetch and process channel.
process process_ensembl_metadata {
  tag "$name ${mysql.db_name}"
  maxForks params.import_metadata.ensembl.max_forks

  input:
  set val(mysql), val(name), file(sql) from ensembl_metadata_tasks

  output:
  file("${name}.csv") optional true into ensembl_output
  val("$name,$mysql.db_name") into ensembl_imported

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
  .mix(term_info)
  .mix(references_output)
  .mix(ensembl_output)
  .mix(ensembl_imported_output)
  .map { f ->
    name = f.getBaseName()
    ctl = file("files/import-data/load/${name.replace('_', '-')}.ctl")
    [[name, ctl], f]
  }
  .filter {
    status = it[0][1].exists()
    if (!status) {
      log.info "Skipping data ${it[0][1].getBaseName()}"
    }
    status
  }
  .groupTuple()
  .map { it -> [it[0][0], it[0][1], it[1]] }
  .set { to_load }

process merge_and_import {
  echo true
  tag name

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
  .into { flag_for_qa; flag_for_mapping; flag_for_feedback }

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
  tag name
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
  .set { flag_for_precompute }
