#!/usr/bin/env nextflow

// This channel will store all the outputs for all data to import. We then split
// the channel into several different channels at the end. Doing this prevents
// us from having many large duplicated blocks. Also adding a new file type to
// import is much easier this way, as we only have to modify in a few places
// instead of every possible process that could produce it.
raw_output = Channel.empty()
to_fetch = Channel.empty()

for (entry in params.import_data.databases) {
  if (!entry.value) {
    continue
  }

  name = entry.key
  database = params.databases[name]

  // The custom property should be a list of directories that contain csv files to
  // import. For example if custom was ./custom we would expect it to contain:
  // ./custom/ac_info.csv, ./custom/genomic_locations.csv, etc. The files may be
  // nested so: ./custom/mirbase/ac_info.csv would also work.
  if (name == "custom") {
    raw_output.mix(Channel.fromPath("${base}/**.csv")).flatten().set { raw_output }
  } else
    db_channel = Channel.value([name, database])
    to_fetch.mix(db_channel).set { to_fetch }
  }
}

//=============================================================================
// Fetch data as well as all extra data for import
//=============================================================================

process fetch_data {
  input:
  set val(name), val(database) from to_fetch

  output:
  set val(name), file("${database.pattern}") into fetched

  script:
  clean_cmd = ''
  for (rm_pattern in database.get('excluded_patterns', [])) {
    clean_cmd = find_cmd + "find . -name '$rm_pattern}' | xargs rm\n"
  }

  fetch_cmd = database.get('fetch_cmd', 'fetch')
  remote = database.get('fetch', '')

  """
  $fetch_cmd '$remote' '${database.pattern}'
  $clean_cmd
  """
}

process fetch_pdb_extra {
  when:
  params.import_data.databases['pdb']

  output:
  set val('pdb'), file('pdb-extra.json') into pdb_extra

  """
  rnac pdb extra pdb-extra.json
  """
}

process fetch_rfam_extra {
  when:
  params.import_data.databases['rfam'] || params.import_data.databases['ensembl']

  input:
  file(query) from Channel.fromPath('files/import-data/rfam/families.sql')

  output:
  set val('rfam'), file('families.tsv') into rfam_extra
  set val('ensembl'), file('families.tsv') into ensembl_extra

  script:
  """
  mysql \
    --host ${params.databases.rfam.mysql.host} \
    --port ${params.databases.rfam.mysql.port} \
    --user ${params.databases.rfam.mysql.user} \
    --database ${params.databases.rfam.mysql.db_name} \
    < $query > families.tsv
  """
}

Channel
  .from(params.databases.ena.tpa_urls)
  .collectFile(name: "urls", newLine: true)
  .set { tpa_url_file }

process fetch_ena_extra {
  when:
  params.import_data.databases['ena']

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
  params.import_data.databases['rgd']

  output:
  set val('rgd'), file('genes.txt') into rgd_extra

  """
  wget -O genes.txt ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt
  """
}

Channel.empty()
  .mix(ena_extra)
  .mix(ensembl_extra)
  .mix(rfam_extra)
  .mix(rgd_extra)
  .mix(pdb_extra)
  .groupBy()
  .set { extra }

fetched
  .flatMap { name, filenames ->
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

//=============================================================================
// Process data
//=============================================================================

process process_data {
  memory params.databases[name].get('memory', '2 GB')

  input:
  set val(name), file(input_file), val(extra) from to_process

  output:
  file "*.csv" into processed_output mode flatten

  script:
  compressed = input_file.toString().endsWith('.gz')
  filename = compressed ? "-" : "$input_file"
  prefix = compressed ? "zcat ${input_file} |" : ''
  """
  $prefix rnac external ${name} ${filename} ${extra_data}
  """
}

//=============================================================================
// Import and release data
//=============================================================================

raw_output
  .mix(processed_output)
  .map { f ->
    filename = f.getName()
    name = filename.take(filename.lastIndexOf('.'))
    ctl = file("files/import-data/${name.replace('_', '-')}.ctl")
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
  set -o pipefail

  mkdir merged
  find . -name 'raw*.csv' |\
  xargs cat |\
  split --additional-suffix=.csv -dC ${params.import_data.chunk_size} - merged/${name}

  cp $ctl merged/$ctl
  cd merged
  pgloader --on-error-stop $ctl
  """
}

loaded
  .map { name -> file("files/import-data/post-release/${name.replace('_', '-')}.sql") }
  .mix(Channel.fromPath('files/import-data/post-release/cleanup.sql'))
  .filter { f -> f.exists() }
  .set { post_scripts }

process release {
  echo true
  maxForks 1

  input:
  file(post) from post_scripts.collect()

  """
  set -o pipefail

  rnac run-release
  find . -name '*.sql' -print0 |\
  sort -z |\
  xargs -r0 -I {} psql -v ON_ERROR_STOP=1 -f {} "$PGDATABASE"
  """
}
