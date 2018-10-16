#!/usr/bin/env nextflow

// This channel will store all the outputs for all data to import. We then split
// the channel into several different channels at the end. Doing this prevents
// us from having many large duplicated blocks. Also adding a new file type to
// import is much easier this way, as we only have to modify in a few places
// instead of every possible process that could produce it.
raw_output = Channel.empty()

// Because we use very generic globs to get files to import we have to filter
// the filenames ot known ones. This maps from the allowed filenames to the
// control files that will be used to load the data.
IMPORTABLE = [
  "accessions.csv": 'files/import-data/accessions.ctl',
  "genomic_locations.csv": 'files/import-data/locations.ctl',
  "seq_long.csv": 'files/import-data/long-sequences.ctl',
  "seq_short.csv": 'files/import-data/short-sequences.ctl',
  "refs.csv": 'files/import-data/references.ctl',
  "secondary_structure.csv": 'files/import-data/secondary.ctl',
  "related_sequences.csv": 'files/import-data/related-sequences.ctl',
  "features.csv": 'files/import-data/features.ctl',
  "sequence_regions.csv": 'files/import-data/regions.ctl',
]

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
  .set { extra }

fetched
  .flatMap { name, filenames ->
    results = []
    filenames.each { filename ->
      results << [name, filename]
    }
    results
  }
  .join(extra, remainder: true)
  .filter { it[0] != null }
  .map { it -> [it[0][0], it[0][1], it[1]] }
  .set { to_import }

//=============================================================================
// Process data
//=============================================================================

process process_data {
  memory params.databases[name].get('memory', '2 GB')

  input:
  set val(name), file(input_file), val(extra) from to_import

  output:
  file "*.csv" into processed_output mode flatten

  script:
  if (input_file.toString().endsWith(".gz")) {
    """
    zcat ${input_file} | rnac external ${name} - ${extra == null ? '' : file(extra)}
    """
  } else {
    """
    rnac external ${name} ${input_file} ${extra == null ? '' : file(extra)}
    """
  }
}

//=============================================================================
// Import and release data
//=============================================================================

raw_output
  .mix(dataless_output)
  .mix(processed_output)
  .filter { f -> IMPORTABLE.keySet().findIndexOf { p -> f.getName() ==~ p } > -1 }
  .map { f ->
    filename = f.getName()
    name = filename.take(filename.lastIndexOf('.'))
    ctl = file(IMPORTABLE[filename])
    [[name, ctl], f]
  }
  .groupTuple()
  .map { it -> [it[0][0], it[0][1], it[1]] }
  .set { to_load }

process merge_and_import {
  echo true

  input:
  set val(name), file(ctl), file('raw*.csv') from to_load

  output:
  file(ctl) into loaded

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

process release {
  echo true
  maxForks 1

  input:
  file('*.ctl') from loaded.collect()
  file(post) from Channel.fromPath('files/import-data/release/post/*.sql').collect()

  """
  rnac run-release
  find . -name '*.sql' -print0 | sort -z | xargs -r0 cat > post-command
  psql -f post-command "$PGDATABASE"
  """
}
