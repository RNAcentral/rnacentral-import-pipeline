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
  if (!entry.value) {
    continue
  }

  if (entry.key == "custom") {
    processed_data.addAll(files("${entry.value}/**.csv"))
    continue
  }

  def db_name = entry.key
  if (!params.databases[db_name]) {
    error "Database ${db_name} is not configured"
  }

  def database = params.databases[db_name]
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

  def linked_to = info.get('linked_databases', [])
  def should_run = info.get('run', false) ||
        (linked_to == '*' && params.import_data.databases) ||
        linked_to.inject(false) { acc, n -> acc || params.import_data.databases[n] }

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

all_processed_output
  .mix(all_fetched_and_processed_output)
  .filter { f -> !f.isEmpty() }
  .choice(processed_output, terms, refs) { f ->
    def names = ["terms.csv", "ref_ids.csv"]
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
  rnac ols lookup-terms unique-terms.txt ontology_terms.csv
  """
}

refs
  .collect()
  .set { refs_to_split }

process merge_and_split_all_publications {
  input:
  file("ref_ids*.csv") from refs_to_split

  output:
  file('split-refs/*.csv') into split_references

  """
  set -o pipefail

  mkdir split-refs
  find . -name 'ref_ids*.csv' | xargs cat | sort -u > all-ids
  split --additional-suffix=".csv" --number l/${params.lookup_publications.maxForks} all-ids split-refs/refs
  """
}

process fetch_publications {
  when { params.import_data.length > 0 }

  output:
  file('out') into refs_database

  """
  fetch europepmc http://europepmc.org/ftp/pmclitemetadata/PMCLiteMetadata.tgz out
  """
}

split_references
  .flatten()
  .combine(refs_database)
  .set { refs_to_lookup }

process lookup_publications {
  memory 4.GB
  maxForks params.lookup_publications.maxForks

  input:
  set file(refs), file(pubs) from refs_to_lookup

  output:
  file("references.csv") into references_output

  script:
  """
  rnac europepmc stream-lookup --ignore-missing --allow-fallback $pubs $refs references.csv
  """
}

process create_load_tables {
  input:
  file(create) from Channel.fromPath('files/schema/create_load.sql')

  output:
  val('done') into created_tables

  """
  psql -v ON_ERROR_STOP=1 -f $create "$PGDATABASE"
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
  .filter { f -> !f.isEmpty() }
  .map { f ->
    def name = f.getBaseName()
    def ctl = file("files/import-data/load/${name.replace('_', '-')}.ctl")
    [[name, ctl], f]
  }
  .filter {
    def status = it[0][1].exists()
    if (!status) {
      log.info "Skipping data ${it[0][1].getBaseName()}"
    }
    status
  }
  .groupTuple()
  .map { it -> [it[0][0], it[0][1], it[1]] }
  .combine(created_tables)
  .set { to_load }

process merge_and_import {
  memory 4.GB
  tag { name }

  input:
  set val(name), file(ctl), file('raw*.csv'), val(flag) from to_load

  output:
  val(name) into (pre_loaded, post_loaded)

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  """
}

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
  maxForks 1

  when:
  data_to_fetch_and_process || data_to_process

  input:
  file(pre_sql) from pre_scripts
  file(post_sql) from post_scripts
  file(limits) from Channel.fromPath('files/import-data/limits.json')

  output:
  val('done') into post_release

  script:
  def should_release = Utils.must_release(params.import_data, params.databases)
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

  ${should_release ? '' : '# ' }rnac check-release $limits
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

process generate_qa_scan_files {
  input:
  set val(name), val(base) from files_to_prepare

  output:
  set val(name), file(name) into qa_scan_files
  set val(name), file('version_file') into qa_version_files

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
    fetch generic "$base/README" version_file
    """
  } else {
    error("Unknown QA to prepare: $name")
  }
}

flag_for_qa
  .combine(Channel.fromPath('files/qa/*.sql').flatten())
  .map { flag, fn -> [fn.getBaseName(), fn] }
  .filter { n, fn -> params.qa[n].run }
  .join(qa_version_files)
  .set { qa_queries }

process fetch_qa_sequences {
  memory 20.GB

  input:
  set val(name), file(query), file(version) from qa_queries

  output:
  set val(name), file('parts/*.fasta') into split_qa_sequences
  set val(name), file('attempted.csv') into qa_track_attempted

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${params.qa[name].chunk_size} --out-dir 'parts/' shuffled.fasta

  rnac qa create-attempted raw.json $name $version attempted.csv
  """
}

split_qa_sequences
  .join(qa_scan_files)
  .flatMap { name, fns, extra -> fns.inject([]) { acc, fn -> acc << [name, fn, extra] } }
  .set { sequences_to_scan }

process qa_scan {
  tag { name }
  cpus { params.qa[name].cpus }
  memory { params.qa[name].memory * params.qa[name].cpus }

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
  .join(qa_track_attempted)
  .map { n, files, ctl, attempted -> 
    [n, files, ctl, attempted, file("files/qa/attempted/${n}.ctl")]
  }
  .set { hits_to_import }

process import_qa_data {
  tag { "qa-$name" }

  input:
  set val(name), file('raw*.csv'), file(ctl), file('attempted*.csv'), file(attempted_ctl) from hits_to_import

  output:
  val("$name done") into qa_imported

  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} $name
  split-and-load $attempted_ctl 'attempted*.csv' ${params.import_data.chunk_size} attempted-$name
  """
}

qa_imported
  .collect()
  .ifEmpty('no qa')
  .set { qa_flag_for_precompute }

//=============================================================================
// Genome mapping
//=============================================================================

process genome_mapping_setup {
  when:
  params.genome_mapping.run

  input:
  val(flag) from flag_for_mapping
  file(setup) from Channel.fromPath('files/genome-mapping/setup.sql')
  file(query) from Channel.fromPath('files/genome-mapping/find-species.sql')

  output:
  file('species.csv') into raw_genomes

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=${params.genome_mapping.to_map_table} \
    -v species_to_map=${params.genome_mapping.species_table} \
    -v min_length=${params.genome_mapping.min_length} \
    -v max_length=${params.genome_mapping.max_length} \
    -f "$setup" "$PGDATABASE"
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=${params.genome_mapping.to_map_table} \
    -v species_to_map=${params.genome_mapping.species_table} \
    -f "$query" "$PGDATABASE" > species.csv
  """
}

raw_genomes
  .splitCsv()
  .filter { s, a, t, d -> !params.genome_mapping.species_excluded_from_mapping.contains(s) }
  .into { assemblies; genomes_to_fetch }

assemblies
  .combine(Channel.fromPath('files/genome-mapping/find-unmapped.sql'))
  .set { assemblies_to_fetch }

process fetch_unmapped_sequences {
  tag { species }
  maxForks params.genome_mapping.fetch_unmapped_sequences.directives.maxForks
  clusterOptions '-sp 100'

  input:
  set val(species), val(assembly_id), val(taxid), val(division), file(query) from assemblies_to_fetch

  output:
  set species, file('parts/*.fasta') into split_mappable_sequences
  file('attempted.csv') into genome_mapping_attempted_sequences

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v taxid=$taxid \
    -v assembly_id=$assembly_id \
    -v tablename=${params.genome_mapping.to_map_table} \
    -v species_to_map=${params.genome_mapping.species_table} \
    -f "$query" \
    "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  split-sequences \
    --max-nucleotides ${params.genome_mapping.fetch_unmapped_sequences.nucleotides_per_chunk} \
    --max-sequences ${params.genome_mapping.fetch_unmapped_sequences.sequences_per_chunk} \
      shuffled.fasta parts

  rnac genome-mapping create-attempted raw.json $assembly_id attempted.csv
  """
}

process download_genome {
  tag { species }
  memory { params.genome_mapping.download_genome.directives.memory }
  errorStrategy 'ignore'

  input:
  set val(species), val(assembly), val(taxid), val(division) from genomes_to_fetch

  output:
  set val(species), val(assembly), file('parts/*.{2bit,ooc}') into genomes

  """
  set -o pipefail

  rnac genome-mapping url-for --host=$division $species $assembly - |\
    xargs -I {} fetch generic '{}' ${species}.fasta.gz

  gzip -d ${species}.fasta.gz
  split-sequences \
    --max-nucleotides ${params.genome_mapping.download_genome.nucleotides_per_chunk} \
    --max-sequences ${params.genome_mapping.download_genome.sequences_per_chunk} \
    ${species}.fasta parts

  find parts -name '*.fasta' |\
    xargs -I {} faToTwoBit -noMask {} {}.2bit

  find parts -name '*.fasta' |\
  xargs -I {} \
    blat \
      -makeOoc={}.ooc \
      -stepSize=${params.genome_mapping.blat.options.step_size} \
      -repMatch=${params.genome_mapping.blat.options.rep_match} \
      -minScore=${params.genome_mapping.blat.options.min_score} \
      {} /dev/null /dev/null
  """
}

genomes
  .join(split_mappable_sequences)
  .flatMap { species, assembly, genome_chunks, chunks ->
    [genome_chunks.collate(2), chunks]
      .combinations()
      .inject([]) { acc, it -> acc << [species, assembly] + it.flatten() }
  }
  .set { targets }

process blat {
  tag { "${species}-${genome.baseName}-${chunk.baseName}" }
  memory { params.genome_mapping.blat.directives.memory }
  errorStrategy 'finish'

  input:
  set val(species), val(assembly), file(genome), file(ooc), file(chunk) from targets

  output:
  set val(species), file('selected.json') into blat_results

  """
  set -o pipefail

  blat \
    -ooc=$ooc \
    -noHead \
    -q=rna \
    -stepSize=${params.genome_mapping.blat.options.step_size} \
    -repMatch=${params.genome_mapping.blat.options.rep_match} \
    -minScore=${params.genome_mapping.blat.options.min_score} \
    -minIdentity=${params.genome_mapping.blat.options.min_identity} \
    $genome $chunk output.psl

  sort -k 10 output.psl |\
    rnac genome-mapping blat serialize $assembly - - |\
    rnac genome-mapping blat select - selected.json
  """
}

 blat_results
  .groupTuple()
  .set { species_results }

process select_mapped_locations {
  tag { species }
  memory { '20GB' }

  input:
  set val(species), file('selected*.json') from species_results

  output:
  file('locations.csv') into selected_locations

  """
  set -o pipefail

  find . -name 'selected*.json' |\
    xargs cat |\
    rnac genome-mapping blat select --sort - - |\
    rnac genome-mapping blat as-importable - locations.csv
  """
}

selected_locations
  .collect()
  .set { blat_to_import }

genome_mapping_attempted_sequences
  .collect()
  .set { genome_mapping_attempted }

process load_genome_mapping {
  maxForks 1

  input:
  file('raw*.csv') from blat_to_import
  file(ctl) from Channel.fromPath('files/genome-mapping/load.ctl')
  file('attempted*.csv') from genome_mapping_attempted
  file(attempted_ctl) from Channel.fromPath('files/genome-mapping/attempted.ctl')

  output:
  val('done') into genome_mapping_status

  script:
  """
  split-and-load $ctl 'raw*.csv' ${params.import_data.chunk_size} genome-mapping
  split-and-load $attempted_ctl 'attempted*.csv' ${params.import_data.chunk_size} genome-mapping-attempted
  psql -v ON_ERROR_STOP=1 -v tablename=params.genome_mapping.to_map_table -c 'DROP TABLE :tablename' $PGDATABASE
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
// Run traveler structures
//=============================================================================

flag_for_secondary
  .combine(Channel.fromPath("files/traveler/find-sequences.sql"))
  .map { flag, query -> query }
  .set { traveler_setup }

process find_possible_traveler_sequences {
  tag { "${model}" }
  memory params.secondary.find_possible.memory
  maxForks params.secondary.find_possible.maxForks
  clusterOptions '-sp 100'

  input:
  file(query) from traveler_setup

  output:
  set file('parts/*.fasta') into to_layout mode flatten

  script:
  def chunk_size = params.secondary.sequence_chunk_size
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v 'model=$model' \
    -v 'tablename=${params.secondary.tablename}' \
    -f "$query" "$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral.fasta
  seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
  seqkit split --two-pass --by-size ${chunk_size} --out-dir 'parts/' shuffled.fasta
  """
}

process layout_sequences {
  tag { "${sequences}" }
  memory params.secondary.layout.memory
  container 'rnacentral/auto-traveler:dev'

  input:
  file(sequences) from to_layout

  output:
  set file("$sequences"), file('output') into secondary_to_parse
  file('output/hits.txt') optional true into secondary_hits
  file('output/*.stk') optional true into secondary_stk

  """
  auto-traveler.py draw $sequences output/
  """
}

process parse_layout {
  input:
  set file('sequences.fasta'), file(to_parse) from secondary_to_parse

  output:
  file("data.csv") into secondary_to_import
  file('attempted.csv') into traveler_attempted_sequences

  """
  rnac traveler process-svgs $family $to_parse data.csv
  rnac traveler create-attempted $sequences attempted.csv
  """
}

secondary_to_import
  .collect()
  .combine(Channel.fromPath("files/traveler/load.ctl"))
  .set { secondary_to_import }

traveler_attempted_sequences
  .collect()
  .combine(Channel.fromPath('files/traveler/attempted.ctl'))
  .set { traveler_attempted }

process store_secondary_structures {
  memory params.secondary.store.memory

  input:
  set file('data*.csv'), file(ctl), file(rna_types_sql) from secondary_to_import
  set file('attempted*.csv'), file(attempted_ctl) from traveler_attempted

  """
  split-and-load $ctl 'data*.csv' ${params.secondary.data_chunk_size} traveler-data
  split-and-load $attempted_ctl 'attemped*.csv' ${params.secondary.data_chunk_size} traveler-attempted
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
