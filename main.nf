#!/usr/bin/env nextflow

assert params.precompute.tablename != 'rna' : "Should not use 'rna' table for precompute"

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

// Create all metadata tasks. We have to run metadata import
for (entry in params.metadata) {
  def name = entry.key
  def info = entry.value

  if (!data_to_fetch && !data_to_fetch_and_process) {
    continue
  }

  def linked_to = info.get('linked_databases', [])
  def should_run = linked_to
    .inject(info.get('always_run', false)) { acc, n -> acc || params.import_data.databases[n] }

  if (!should_run) {
    continue
  }

  def source = DataSource.build("metadata-$name", info)
  source.process.group = "/rnc/metadata/$name"
  if (DataSource.is_parallel_task(source)) {
    data_to_fetch.addAll(source.inputs)
    data_to_process << [source.name, source]
  } else {
    data_to_fetch_and_process << source
  }
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

//=============================================================================
// Import and release data
//=============================================================================

raw_output
  .mix(
    processed_output,
    term_info,
    references_output,
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
  .toList()
  .set { pre_scripts }

post_loaded
  .flatMap { n -> file("files/import-data/post-release/*__${n.replace('_', '-')}.sql") }
  .mix(
    Channel.fromPath([
      'files/import-data/post-release/000__populate_precompute.sql',
      'files/import-data/post-release/999__cleanup.sql',
    ])
  )
  .filter { f -> f.exists() }
  .toList()
  .set { post_scripts }

process release {
  echo true
  maxForks 1

  when:
  data_to_fetch_and_process || data_to_process

  input:
  file(pre_sql) from pre_scripts
  file(post_sql) from post_scripts

  output:
  val('done') into post_release

  script:
  def should_release = params.import_data.databases.inject(false) { s, e -> s || e.value }
  def pre = file("work/pre-release")
  def post = file("work/post-release")
  """
  set -o pipefail

  run_sql() {
    local fn="\$1"
    if [[ -s "\$fn" ]]; then
      while IFS='' read -r "script" || [[ -n "\$script" ]]; do
        if [[ ! -z "\$script" ]]; then
          echo "Running: \$fn/\$script"
          psql -v ON_ERROR_STOP=1 -f \$script "$PGDATABASE"
        fi
      done < "\$fn"
    fi
  }

  run_sql "${ Utils.write_ordered(pre, pre_sql.inject([]) { a, fn -> a << fn.getName() }) }"
  ${should_release ? '' : '# ' }rnac run-release
  run_sql "${ Utils.write_ordered(post, post_sql.inject([]) { a, fn -> a << fn.getName() }) }"
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
  .map { flag, fn -> [flag, fn.getBaseName(), fn] }
  .filter { f, n, fn -> params.qa[n].run }
  .set { qa_queries }

process fetch_qa_sequences {

  input:
  set val(status), val(name), file(query) from qa_queries

  output:
  set val(name), file('parts/*.fasta') into split_qa_sequences

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${params.qa[name].chunk_size} --out-dir 'parts/' shuffled.fasta
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
    fetch generic "$base/Rfam.cm.gz" Rfam.cm.gz
    gzip -d *.gz
    cmpress Rfam.cm
    cd ..
    """
  } else {
    error("Unknown QA to prepare: $name")
  }
}

split_qa_sequences
  .join(qa_scan_files)
  .flatMap { name, fns, extra -> fns.inject([]) { acc, fn -> acc << [name, fn, extra] } }
  .set { sequences_to_scan }

process qa_scan {
  tag { name }
  cpus { params.qa[name].cpus }
  memory { params.qa[name].memory }

  input:
  set val(name), file('sequences.fasta'), file(dir) from sequences_to_scan

  output:
  set val(name), file('hits.csv') into qa_scan_results

  script:
  if (name == 'rfam') {
    """
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
  .map { n, files -> [n, files, file("files/qa/${n}.ctl")] }
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
  val("$name done") into qa_imported

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  """
}

qa_imported
  .collect()
  .ifEmpty('no qa')
  .set { qa_flag_for_precompute }

//=============================================================================
// Genome mapping
//=============================================================================

process species_to_map {
  executor 'local'

  when:
  params.genome_mapping.run

  input:
  val(flag) from flag_for_mapping
  file(query) from Channel.fromPath('files/genome-mapping/mappable.sql')

  output:
  stdout into raw_genomes

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE"
  """
}

raw_genomes
  .splitCsv()
  .filter { s, a, t, d -> !params.genome_mapping.species_excluded_from_mapping.contains(s) }
  .into { assemblies; genomes_to_fetch; assembly_tracking }

assemblies
  .combine(Channel.fromPath('files/genome-mapping/find-unmapped.sql'))
  .set { assemblies_to_fetch }

process fetch_unmapped_sequences {
  tag { species }
  scratch true
  maxForks 5
  errorStrategy 'ignore'

  input:
  set val(species), val(assembly_id), val(taxid), val(division), file(query) from assemblies_to_fetch

  output:
  set species, file('parts/*.fasta') into split_mappable_sequences

  script:
  """
  psql -v ON_ERROR_STOP=1 -v taxid=$taxid -v assembly_id=$assembly_id -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size "${params.genome_mapping.chunk_size}" --out-dir 'parts/' shuffled.fasta
  """
}

process download_genome {
  tag { species }
  memory 30.GB
  scratch true

  input:
  set val(species), val(assembly), val(taxid), val(division) from genomes_to_fetch

  output:
  set val(species), file('*.fa'), file('11.ooc') into genomes

  script:
  def engine = new groovy.text.SimpleTemplateEngine()
  def url = engine
    .createTemplate(params.genome_mapping.sources[division])
    .make([species: species])
    .toString()
  """
  fetch genome '$url' ${species}.fasta

  blat \
    -makeOoc=11.ooc \
    -stepSize=${params.genome_mapping.blat_options.step_size} \
    -repMatch=${params.genome_mapping.blat_options.rep_match} \
    -minScore=${params.genome_mapping.blat_options.min_score} \
    ${species}.fasta /dev/null /dev/null
  """
}

genomes
  .join(split_mappable_sequences)
  .flatMap { species, chrs, ooc_file, chunks ->
    [chrs, chunks].combinations().inject([]) { acc, it -> acc << [species, ooc_file] + it }
  }
  .filter { s, o, c, t -> !c.empty() }
  .filter { s, o, c, t -> !params.genome_mapping.chromosomes_excluded_from_mapping.contains(c.getBaseName()) }
  .set { targets }

process blat {
  memory 10.GB
  errorStrategy 'finish'

  input:
  set val(species), file(ooc), file(chromosome), file(chunk) from targets

  output:
  set val(species), file('output.psl') into blat_results

  """
  blat \
    -ooc=$ooc \
    -noHead \
    -q=rna \
    -stepSize=${params.genome_mapping.blat_options.step_size} \
    -repMatch=${params.genome_mapping.blat_options.rep_match} \
    -minScore=${params.genome_mapping.blat_options.min_score} \
    -minIdentity=${params.genome_mapping.blat_options.min_identity} \
    $chromosome $chunk output.psl
  """
}

 blat_results
  .groupTuple()
  .join(assembly_tracking)
  .map { species, psl, assembly_id, taxid, division -> [psl, species, assembly_id] }
  .set { species_results }

process select_mapped_locations {
  tag { species }
  memory '15 GB'

  input:
  set file('output*.psl'), val(species), val(assembly_id) from species_results

  output:
  file 'locations.csv' into selected_locations

  """
  set -o pipefail

  sort -k 10 output*.psl > sorted.psl
  rnac genome-mapping select-hits $assembly_id sorted.psl locations.csv
  """
}

process load_genome_mapping {
  maxForks 1

  input:
  file('raw*.csv') from selected_locations.collect()
  file(ctl) from Channel.fromPath('files/genome-mapping/load.ctl')

  output:
  val('done') into genome_mapping_status

  script:
  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} genome-mapping
  """
}

genome_mapping_status
  .ifEmpty('no mapping')
  .set { mapping_flag_for_precompute }

//=============================================================================
// Run precompute of selected data
//=============================================================================

qa_flag_for_precompute
  .combine(mapping_flag_for_precompute)
  .map { qa, gm -> file("files/precompute/methods/${params.precompute.method.replace('_', '-')}.sql") }
  .set { precompute_upi_queries }

process find_precompute_upis {
  when: params.precompute.run

  input:
  file(sql) from precompute_upi_queries

  output:
  file('ranges.txt') into raw_ranges

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=$params.precompute.tablename \
    -f "$sql" "$PGDATABASE"
  rnac upi-ranges --table-name $params.precompute.tablename ${params.precompute.max_entries} ranges.txt
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

  output:
  val('precompute') into post_precompute

  script:
  def tablename = params.precompute.tablename
  """
  split-and-load $pre_ctl 'precompute*.csv' ${params.import_data.chunk_size} precompute
  split-and-load $qa_ctl 'qa*.csv' ${params.import_data.chunk_size} qa
  psql -v ON_ERROR_STOP=1 -f $post "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -c 'DROP TABLE IF EXISTS $tablename' "$PGDATABASE"
  """
}

post_precompute
  .ifEmpty('no precompute')
  .into { flag_for_feedback; flag_for_secondary }

//=============================================================================
// Compute secondary structures
//=============================================================================

flag_for_secondary
  .combine(Channel.fromPath("files/secondary-structures/find-sequences.sql"))
  .set { secondary_query }

process find_possible_secondary_sequences {
  when:
  params.secondary.run

  input:
  set val(flag), file(query) from secondary_query

  output:
  file('parts/*.fasta') into sequences_to_ribotype mode flatten

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${params.secondary.sequence_chunk_size} --out-dir 'parts/' shuffled.fasta
  """
}

process fetch_traveler_data {
  when:
  params.secondary.run

  output:
  set file('auto-traveler/data/cms/'), file('auto-traveler/data/crw-fasta-no-pseudoknots/'), file('auto-traveler/data/crw-ps/') into traveler_data

  """
  git clone https://github.com/RNAcentral/auto-traveler.git
  cd auto-traveler
  git checkout "${params.secondary.auto_traveler_version}"
  wget -O cms.tar.gz '${params.secondary.cm_library}'
  # We are going to ignore some errors due to using mac tar to build the tarball
  tar xf cms.tar.gz
  python utils/generate_model_info.py --cm-library data/cms
  """
}

sequences_to_ribotype
  .combine(traveler_data)
  .set { to_layout }

process layout_sequences {
  input:
  set file(sequences), file(cm), file(fasta), file(ps) from to_layout

  output:
  file("output/") into secondary_to_import

  """
  auto-traveler.py --cm-library $cm --fasta-library $fasta --ps-library $ps $sequences output/
  """
}

secondary_to_import
  .collect()
  .combine(Channel.fromPath("files/secondary-structures/load.ctl"))
  .set { secondary_to_import }

process store_secondary_structures {
  input:
  set file(svg_dir), file(ctl) from secondary_to_import

  """
  rnac secondary process-svgs $svg_dir data.csv
  split-and-load $ctl data.csv ${params.secondary.data_chunk_size} traveler-data
  """
}

//=============================================================================
// Compute feedback reports
//=============================================================================

flag_for_feedback
  .combine(Channel.fromPath('files/precompute/find-mod-info.sql'))
  .set { feedback_queries }

process mods_for_feedback {
  when:
  params.feedback.run

  input:
  set val(status), file(query) from feedback_queries

  output:
  file('info') into raw_mods

  script:
  names = []
  for (name in params.feedback.databases) {
    names << "'${name.toUpperCase()}'"
  }
  names = '(' + names.join(', ') + ')'
  """
  psql -v ON_ERROR_STOP=1 -v "names=${names}" -f "$query" "$PGDATABASE" > info
  """
}

raw_mods
  .splitCsv()
  .combine(Channel.fromPath('files/ftp-export/genome_coordinates/query.sql'))
  .set { mods }

process generate_feedback_report {
  tag { "${mod}-${assembly}" }
  memory { params.feedback.report.memory * task.attempt }
  errorStrategy 'retry'
  maxRetries 4

  input:
  set val(assembly), val(mod), file(query) from mods

  output:
  file('combined.tsv') into feedback

  """
  find-overlaps $query complete-${mod}.bed $assembly
  """
}

process import_feedback {
  echo true

  input:
  file('feedback*.tsv') from feedback.collect()
  file(ctl) from Channel.fromPath('files/precompute/feedback.ctl')

  """
  split-and-load $ctl 'feedback*.tsv' ${params.import_data.chunk_size} feedback tsv
  """
}

workflow.onComplete {
  if (params.notify) {
    def success = (workflow.success ? '--success' : '--failure');
    def summary = "${workflow.scriptName} completed ${workflow.success ? 'successfully' : 'with errors'} at ${workflow.complete}";
    def msg_file = File.createTempFile("msg", ".txt");
    msg_file << workflow.errorReport ?: 'No errors';
    def cmd = ['bin/slack', 'pipeline-done', '--url', params.notify_url, success, "'$summary'", msg_file.getAbsolutePath()];
    def process = cmd.execute();
    process.waitForProcessOutput(System.out, System.err);
  }
}
