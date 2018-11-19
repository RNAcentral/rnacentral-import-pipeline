#!/usr/bin/env nextflow

def any_database(String... names) { names.any { n -> params.import_data.databases[n] } }

def as_mysql_cmd = { db ->
  def rest = db.db_name ? " --database ${db.db_name}" : ''
  "mysql --host ${db.host} --port ${db.port} --user ${db.user} $rest"
}

// ===========================================================================
// Compute initial tasks
// ===========================================================================

processed_data = []

// This contains what remote data sources we have to fetch. It focuses on the
// main data files, not the 'extra' ones. Separating the two makes rejoining
// them later much easier.
data_to_fetch = []
data_to_process = []

// This contains the tasks for jobs that can be fetched and processed in a
// single task. Some databases, like all JSON based import ones provide a single
// file that we need to process. In that case we don't have to split the
// fetching and processing into separate tasks. We can do 1 task instead of 2
// (makes the job scheduler and ITS happier, and may be faster if the cluser is
// hunder heavy load) and perform a single task that does both.
data_to_fetch_and_process = []

// Here we generate the data that will go into the channels. This is a bit
// cleaner to work with over constantly modifying the channels in the loop.
for (entry in params.import_data.databases) {
  def db_name = entry.key
  def database = params.databases[db_name]

  if (!entry.value || !params.databases[db_name]) {
    continue
  }

  if (db_name == "custom") {
    processed_data.addAll(files("${entry.value}/**.csv"))
    continue
  }

  def source = DataSource.build(db_name, database)
  if (DataSource.is_parallel_task(source)) {
    data_to_fetch.addAll(source.inputs)
    data_to_process << [source.name, source]
  } else {
    data_to_fetch_and_process << source
  }
}

Channel.fromPath('files/import-data/ensembl/*.sql')
  .map { f -> f.getBaseName() }
  .collectFile(name: "possible-data.txt", newLine: true)
  .set { ensembl_possible }

Channel.from(params.metadata.ensembl.mysql, params.metadata.ensembl_genomes.mysql)
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

// If we are processing Rfam or Ensembl we need to update the Rfam metadata.
// These databases use Rfam in one way or another so we should stay up-to-date
// with the current Rfam data. For example they may be creating annotations of a
// partial lncRNA which we do not want to import. By tracking the latest Rfam we
// will know what is a partial lncRNA. This is done by querying the Rfam MySQL
// database and then using the Rfam commands to process the query results. We
// can then import the processed data as normal.
if (any_database('rfam', 'ensembl')) {
  file("files/import-data/rfam/*.sql").each { query ->
    def name = query.getBaseName()
    def input = [produces: "data.tsv", query: query.toString(), command: 'mysql']
    data_to_fetch_and_process << DataSource.build("rfam-$name", [
      inputs: [data: input + params.metadata.rfam.mysql],
      process: [command: "rnac rfam $name"]
    ])
  }
}

// Add the metadata tasks that must always be run when importing data
if  (data_to_fetch || data_to_fetch_and_process) {
  data_to_fetch_and_process << DataSource.build('pub-info', params.metadata.europepmc)
  data_to_fetch_and_process << DataSource.build('ncbi-taxonomy', params.metadata.taxonomy)
}

// Setup the QA files for hmmpress or whatever is needed.
to_prepare = params.qa
  .findAll { key, value -> value.get('run', true) }
  .inject([]) { acc, entry -> acc << [entry.key, entry.value.files] }

// Now setup the channels with all data
Channel.from(processed_data).set { raw_output }
Channel.from(data_to_fetch).set { to_fetch }
Channel.from(data_to_process).set { process_specs }
Channel.from(data_to_fetch_and_process).set { to_fetch_and_process }
Channel.from(to_prepare).set { files_to_prepare }

//=============================================================================
// Fetch data as well as all extra data for import
//=============================================================================

process fetch_data {
  tag { inp.directives.tag }
  memory { inp.directives.memory }
  clusterOptions { "-g $inp.directives.group" }

  input:
  val(inp) from to_fetch

  output:
  set val(inp), file("${inp.produces}") into fetched

  script:
  """
  set -o pipefail

  ${Input.script(inp)}
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

fetched
  .map { t, fs ->
    def ps = t.exclude
    def selected = [fs].flatten().findAll { f -> !Input.is_excluded(t, f.getName()) }
    [t.source, t, selected]
  }
  .groupTuple(by: 0)
  .combine(process_specs, by: 0)
  .flatMap { name, tasks, inputs, spec ->
    [tasks, inputs]
      .transpose()
      .sort { it[0].index }
      .inject([]) { a, d -> a << d[1] }
      .combinations()
      .inject([]) { a, files -> a << [spec, files] }
  }
  .set { to_process }

//=============================================================================
// Process data
//=============================================================================

process fetch_and_process_data {
  tag { spec.process.directives.tag }
  memory { spec.process.directives.memory }
  clusterOptions { "-g $spec.process.directives.group" }

  input:
  val(spec) from to_fetch_and_process

  output:
  file("${spec.process.produces}") into all_fetched_and_processed_output mode flatten

  script:
  def input_files = spec.inputs.inject([]) { acc, inp -> acc << inp.produces }
  def fetches = spec.inputs.inject([]) { acc, inp -> acc << Input.script(inp) }
  """
  set -o pipefail

  ${fetches.join('\n')}
  ${DataSource.process_script(spec, input_files)}
  """
}

process process_data {
  tag { spec.process.directives.tag }
  memory { spec.process.directives.memory }
  clusterOptions { "-g $spec.process.directives.group" }

  input:
  set val(spec), file(input_files) from to_process

  output:
  file("${spec.process.produces}") into all_processed_output mode flatten

  script:
  def filenames = input_files.inject([]) { a, fn -> a << fn.getName() }
  """
  set -o pipefail

  ${DataSource.process_script(spec, filenames)}
  """
}

processed_output = Channel.create()
refs = Channel.create()
terms = Channel.create()
refs_database = Channel.create()

all_processed_output
  .mix(all_fetched_and_processed_output)
  .choice(processed_output, terms, refs, refs_database) { f ->
    names = ["terms.csv", "ref_ids.csv", params.metadata.europepmc.produces]
    return names.indexOf(f.getName()) + 1
  }

process batch_lookup_ontology_information {
  input:
  file('terms*.csv') from terms.collect()

  output:
  file('ontology_terms.csv') into term_info

  script:
  """
  set -o pipefail

  find . -name 'terms*.csv' | xargs cat | sort -u >> unique-terms.txt
  rnac ontologies lookup-terms unique-terms.txt ontology_terms.csv
  """
}

refs.collect().set { refs_to_lookup }

process batch_lookup_publications {
  input:
  file("ref_ids*.csv") from refs_to_lookup
  file("references.db") from refs_database

  output:
  file('references.csv') into references_output

  """
  set -o pipefail

  find . -name 'ref_ids*.csv' | xargs cat | sort -u >> all-ids
  rnac europepmc lookup --allow-fallback references.db all-ids references.csv
  """
}

// This must be it's own task becuase we need to have a specific number of forks
// and this cannot be set dynamically. Otherwise we could just put it into the
// generic fetch and process channel.
process process_ensembl_metadata {
  tag "$name ${mysql.db_name}"
  maxForks params.metadata.ensembl.max_forks
  clusterOptions { "-g /rnc/process/ensembl-metadata" }

  input:
  set val(mysql), val(name), file(sql) from ensembl_metadata_tasks

  output:
  file("${name}.csv") into ensembl_output
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
  .mix(
    processed_output,
    term_info,
    references_output,
    ensembl_output,
    ensembl_imported_output
  )
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
  tag { name }

  input:
  set val(name), file(ctl), file('raw*.csv') from to_load

  output:
  val(name) into loaded

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  """
}

loaded.into { pre_loaded; post_loaded }

pre_loaded
  .flatMap { n -> file("files/import-data/pre-release/*__${n.replace('_', '-')}.sql") }
  .filter { f -> f.exists() }
  .toSortedList()
  .set { pre_scripts }

post_loaded
  .flatMap { n -> file("files/import-data/post-release/*__${n.replace('_', '-')}.sql") }
  .mix(
    Channel.fromPath('files/import-data/post-release/000__populate_precompute.sql'),
    Channel.fromPath('files/import-data/post-release/999__cleanup.sql')
  )
  .filter { f -> f.exists() }
  .toSortedList()
  .set { post_scripts }

process release {
  echo true
  maxForks 1

  input:
  file(pre_sql) from pre_scripts
  file(post_sql) from post_scripts

  output:
  val('done') into post_release

  script:
  """
  set -o pipefail

  run_sql() {
    local fn="\$1"
    while IFS='' read -r "script" || [[ -n "\$script" ]]; do
      echo "Running: \$fn/\$script"
      psql -v ON_ERROR_STOP=1 -f \$script "\$PGDATABASE"
    done < "\$fn"
  }

  echo "${pre_sql.join('\n')}" > pre-release
  echo "${post_sql.join('\n')}" > post-release

  run_sql pre-release
  rnac run-release
  run_sql "post-release"
  """
}

post_release
  .ifEmpty('no release')
  .into { flag_for_qa; flag_for_mapping }

//=============================================================================
// QA scans
//=============================================================================

flag_for_qa
  .combine(Channel.fromPath('files/qa/*.sql').flatten())
  .map { flag, fn ->
    name = fn.getBaseName()
    [flag, name.take(name.indexOf('.')), fn]
  }
  .filter { f, n, fn -> params.qa[n].run }
  .set { qa_queries }

process fetch_sequences {

  input:
  set val(status), val(name), file(query) from qa_queries

  output:
  set val(name), file('parts/*.fasta') into split_sequences mode flatten

  """
  psql -v ON_ERROR_STOP=1 $pg_args -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${params.qa.rfam_scan.chunk_size} --out-dir 'parts/' shuffled.fasta
  """
}

process generate_qa_scan_files {
  input:
  set val(name), val(base) from files_to_prepare

  output:
  set val(name), file(name) into qa_scan_files

  script:
  if (name == "pfam") {
    """
    mkdir $name
    cd $name
    fetch generic "$base/Pfam-A.hmm.gz" Pfam-A.hmm.gz
    fetch generic "$base/Pfam-A.dat.hmm.gz" Pfam-A.dat.hmm.gz
    fetch generic "$base/active_site.dat.gz" active_site.dat.gz
    gzip -d *.gz
    hmmpress Pfam-A.hmm
    cd ..
    """
  } else if (name == "dfam")  {
    """
    mkdir $name
    cd $name
    fetch generic "$base/Dfam.hmm.gz" Dfam.hmm.gz
    gzip -d Dfam.hmm.gz
    hmmpress Dfam.hmm
    cd ..
    """
  } else if (name == "rfam") {
    """
    mkdir $name
    cd $name
    fetch generic "$base/Rfam.clanin" Rfam.clanin
    fetch generic "$base/Rfam.cm.gz Rfam.cm.gz
    gzip -d *.gz
    cmpress Rfam.cm
    cd ..
    """
  } else {
    error("Unknown QA to prepare: $name")
  }
}

split_sequences
  .join(qa_scan_files)
  .flatMap { name, fns, extra -> fns.inject([]) { fn -> [name, fn, extra] } }
  .set { sequences_to_scan }

process qa_scan {
  tag { name }
  cpus { spec.cpus }
  queue { spec.queue }
  module { spec.get('module', '') }
  clusterOptions {
    "-M ${spec.memory} -R 'rusage[mem=${spec.memory}]' ${spec.get('options', '')}"
  }

  input:
  set val(name), file('sequences.fasta'), file(dir) from sequences_to_scan

  output:
  set val(name), file('hits.csv') into qa_scan_results

  script:
  if (name == 'rfam') {
    """
    mpiexec -mca btl ^openbib -np ${params.qa.rfam_scan.cpus} \
    cmscan \
      -o output.inf \
      --tblout results.tblout \
      --clanin $dir/Rfam.clanin \
      --oclan \
      --fmt 2 \
      --acc \
      --cut_ga \
      --rfam \
      --notextw \
      --nohmmonly \
      --mpi \
      "$dir/Rfam.cm" \
      sequences.fasta
    rnac qa $name results.tblout hits.csv
    """
  } else if (name == 'pfam') {
    """
    pfam_scan.pl \
      -fasta sequences.fasta \
      -dir "$dir" \
      -cpus ${params.qa[name].cpus} \
      -outfile raw.tsv
    rnac qa $name raw.tsv hits.csv
    """
  } else if (name == 'dfam') {
    """
    dfamscan.pl \
      -fastafile sequences.fasta \
      --hmmfile "$dir/Dfam.hmm" \
      --cut_ga \
      --cpu ${params.qa[name].cpus} \
      --dfam_outfile raw.txt
    rnac qa $name raw.txt hits.csv
    """
  } else {
    error("Unknown type of QA scan: $name")
  }
}

qa_scan_results
  .groupTuple()
  .map { n, files -> [n, files, Channel.fromPath("files/qa/$n.ctl")] }
  .filter { n, fs, ctl ->
    def status =  ctl.exists()
    if (!status) {
      log.error "Skipping QA results $n"
    }
    status
  }
  .set { hits_to_import }

process import_qa_data {
  tag { "qa-$name" }
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

flag_for_precompute
  .map { flag -> [flag, file("files/precompute/methods/${params.precompute.method.replace('_', '-')}.sql")] }
  .set { precompute_upi_queries }

//=============================================================================
// Run precompute of selected data
//=============================================================================

process find_precompute_upis {
  input:
  set val(flag), file(sql) from precompute_upi_queries

  output:
  file('ranges.txt') into raw_ranges

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$sql" "$PGDATABASE"
  rnac upi-ranges --table-name upis_to_precompute ${params.precompute.max_entries} ranges.txt
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/precompute/query.sql'))
  .set { ranges }

process precompute_range_query {
  tag { "$min-$max" }
  beforeScript 'slack db-work precompute-range || true'
  afterScript 'slack db-done precompute-range || true'
  maxForks params.precompute.maxForks

  input:
  set val(tablename), val(min), val(max), file(query) from ranges

  output:
  set val(min), val(max), file('raw-precompute.json') into precompute_raw

  """
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=${params.precompute.tablename} \
    --variable min=$min \
    --variable max=$max \
    -f "$query" \
    '$PGDATABASE' > raw-precompute.json
  """
}

process precompute_range {
  tag { "$min-$max" }
  memory params.precompute.range.memory

  input:
  set val(min), val(max), file(raw) from precompute_raw

  output:
  file 'precompute.csv' into precompute_results
  file 'qa.csv' into qa_results

  """
  rnac precompute from-file $raw
  """
}

process load_precomputed_data {
  echo true

  beforeScript 'slack db-work loading-precompute || true'
  afterScript 'slack db-done loading-precompute || true'

  input:
  file('precompute*.csv') from precompute_results.collect()
  file('qa*.csv') from qa_results.collect()
  file pre_ctl from Channel.fromPath('files/precompute/load.ctl')
  file qa_ctl from Channel.fromPath('files/precompute/qa.ctl')
  file post from Channel.fromPath('files/precompute/post-load.sql')

  script:
  def tablename = params.precompute.tablename
  """
  split-and-load $pre_ctl 'precompute*.csv' ${params.import_data.chunk_size} precompute
  split-and-load $qa_ctl 'qa*.csv' ${params.import_data.chunk_size} qa
  psql -v ON_ERROR_STOP=1 -f $post "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -c 'DROP TABLE IF EXISTS $tablename' "$PGDATABASE"
  """
}

//=============================================================================
// Compute feedback reports
//=============================================================================

flag_for_feedback
  .combine(Channel.fromPath('files/precompute/find-mod-info.sql'))
  .set { feedback_queries }

// process mods_for_feedback {
//   input:
//   set val(status), file(query) from feedback_queries

//   output:
//   file('info') into raw_mods

//   script:
//   names = []
//   for (name in params.precompute.feedback.databases) {
//     names << "'${name.toUpperCase()}'"
//   }
//   names = '(' + names.join(', ') + ')'
//   """
//   psql -v ON_ERROR_STOP=1 -v "names=${names}" -f "$query" "$PGDATABASE" > info
//   """
// }

// raw_mods
//   .splitCsv()
//   .combine(Channel.fromPath('files/ftp-export/genome_coordinates/query.sql'))
//   .set { mods }

// process generate_feedback_report {
//   memory params.feedback.report.memory

//   input:
//   set val(assembly), val(mod), file(query) from mods

//   output:
//   file('combined.tsv') into feedback

//   """
//   find-overlaps $query complete-${mod}.bed $assembly
//   """
// }

// process import_feedback {
//   echo true

//   input:
//   file('feedback*.tsv') from feedback.collect()
//   file(ctl) from Channel.fromPath('files/precompute/feedback.ctl')

//   """
//   split-and-load $ctl 'feedback*.tsv' ${params.import_data.chunk_size} merged
//   """
// }
