#!/usr/bin/env nextflow

// This channel will store all the outputs for all data to import. We then split
// the channel into several different channels at the end. Doing this prevents
// us from having many large duplicated blocks. Also adding a new file type to
// import is much easier this way, as we only have to modify in a few places
// instead of every possible process that could produce it.
raw_output = Channel.empty()

// Because we use very generic globs to get files to import we have to filter
// the filenames ot known ones. This is the list of known file names. The list
// is also used to shorten the implementation of the .choice method below.
IMPORTABLE_FILENAMES = [
  'ac_info.csv',
  'genomic_locations.csv',
  'seq_long.csv',
  'seq_short.csv',
  'refs.csv',
  'secondary_structure.csv',
  'related_sequences.csv',
  'features.csv',
]

schema_imports = Channel.empty()
dataless_imports = Channel.empty()
remotes_to_fetch = Channel.empty()
locals_to_fetch = Channel.empty()

for (entry in params.import_data.databases) {
  if (!entry.value) {
    continue
  }

  name = entry.key
  database = params.databases[name]

  // The custom property should be a list of directories that contain csv files to
  // import. For example if custom was ./mirbase we would expect it to contain:
  // ./mirbase/ac_info.csv, ./mirbase/genomic_locations.csv, etc. The files may be
  // nested so: ./mirbase/something/ac_info.csv would also work.
  if (name == "custom") {
    raw_output.mix(Channel.fromPath("${base}/**.csv")).flatten().set { raw_output }
  } else if (database.containsKey('json_schema')) {
    schema_imports
      .mix(Channel.from(database.json_schema))
      .set { schema_imports }
  } else if (database.containsKey('remote') && database.remote) {
    db_channel = Channel.value([name, database.remote, database.pattern])
    if (database.remote.contains('http') || database.remote.contains('ftp')) {
      remotes_to_fetch.mix(db_channel).set { remotes_to_fetch }
    } else {
      locals_to_fetch.mix(db_channel).set { locals_to_fetch }
    }
  } else {
    dataless_imports.mix(Channel.from(name)).set { dataless_imports }
  }
}

process json_schema_import {
  memory '4 GB'

  input:
  val(url) from schema_imports

  output:
  file "*.csv" into raw_json_output mode flatten

  script:
  if (url.contains('.gz')) {
    """
    curl "$url" | gzip -d | rnac external json-schema -
    """
  } else {
    """
    curl "$url" | rnac external json-schema -
    """
  }
}

raw_output.mix(raw_json_output).set { raw_output }

process dataless {
  input:
  val(name) from dataless_imports

  output:
  file "*.csv" into raw_dataless_output mode flatten

  script:
  """
  rnac external $name
  """
}

raw_output.mix(raw_dataless_output).set { raw_output }

process remote_fetch {
  input:
  set val(name), val(remote), val(pattern) from remotes_to_fetch

  output:
  set val(name), file("${pattern}") into remote_fetched

  """
  wget "${remote}"
  """
}

process local_fetch {
  input:
  set val(name), val(directory), val(pattern) from locals_to_fetch

  output:
  set val(name), file("${pattern}") into local_fetched

  """
  find ${directory} -name '${pattern}' | xargs -I {} cp {} .
  """
}

metadataless = Channel.create()
with_metadata = Channel.create()
remote_fetched
  .mix(local_fetched)
  .transpose()
  .choice(metadataless, with_metadata) { it ->
    params.databases[it[0]].get('metadata', false) ? 1 : 0
  }

process external_without_metadata {
  input:
  set val(name), file(filename) from metadataless

  output:
  file "*.csv" into raw_metadataless_output mode flatten

  script:
  if (".gz" in filename) {
    """
    zcat $filename | rnac external ${name} -
    """
  } else {
    """
    rnac external ${name} ${filename}
    """
  }
}

raw_output.mix(raw_metadataless_output). set { raw_output }

process fetch_rfam_metadata {
  input:
  file(query) from Channel.fromPath('files/import-data/rfam/families.sql')

  output:
  set val('rfam'), file('families.tsv') into raw_rfam_metadata

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

raw_rfam_metadata.into { rfam_metadata; ensembl_metadata }

ensembl_metadata
  .map { ['ensembl', it[1]] }
  .set { ensembl_metadata }

Channel
  .from(params.databases.ena.tpa_urls)
  .collectFile(name: "urls", newLine: true)
  .set { tpa_url_file }

process fetch_ena_metadata {
  input:
  file tpa_file from tpa_url_file

  output:
  set val('ena'), file("tpa.tsv") into ena_metadata

  """
  cat $tpa_file | xargs wget -O - > tpa.tsv
  """
}

Channel.empty()
  .mix(ena_metadata)
  .mix(ensembl_metadata)
  .mix(rfam_metadata)
  .cross(with_metadata)
  .map { meta, data -> data + meta[1..-1] }
  .set { metadata }

process import_with_metadata {
  input:
  set val(name), file(input_file), file("metadata*") from metadata

  output:
  file "*.csv" into raw_with_metadata_output mode flatten

  script:
  if (input_file.contains(".gz")) {
    """
    zcat ${input_file} | rnac external ${name} - metadata*
    """
  } else {
    """
    rnac external ${name} ${input_file} metadata*
    """
  }
}

raw_output.mix(raw_with_metadata_output).set { raw_output }

accessions = Channel.create()
locations = Channel.create()
long_sequences = Channel.create()
short_sequences = Channel.create()
refs = Channel.create()
secondary = Channel.create()
related = Channel.create()
features = Channel.create()

raw_output
  .filter { filename -> IMPORTABLE_FILENAMES.contains(filename.getName()) }
  .choice(
    accessions,
    locations,
    long_sequences,
    short_sequences,
    refs,
    secondary,
    related,
    features,
  ) { filename -> IMPORTABLE_FILENAMES.indexOf(filename.getName()) }

process pgload_data {

  input:
  file 'ac*' from accessions.collect()
  file 'location*' from locations.collect()
  file 'long*' from long_sequences.collect()
  file 'ref*' from refs.collect()
  file 'short*' from short_sequences.collect()

  file acc_ctl from Channel.fromPath('files/import-data/accessions.ctl')
  file ref_ctl from Channel.fromPath('files/import-data/references.ctl')
  file coord_ctl from Channel.fromPath('files/import-data/coordinates.ctl')
  file short_ctl from Channel.fromPath('files/import-data/short-sequences.ctl')
  file long_ctl from Channel.fromPath('files/import-data/long-sequences.ctl')

  """
  mkdir ac_info genome_locations long refs short secondary_structures

  # cat ac* | split -dC ${params.import_data.chunk_size} - ac_info/chunk_
  # cat location* | split -dC ${params.import_data.chunk_size} - genome_locations/chunk_
  # cat long* | split -dC ${params.import_data.chunk_size} - long/chunk_
  # cat ref* | split -dC ${params.import_data.chunk_size} - refs/chunk_
  # cat short* | split -dC ${params.import_data.chunk_size} - short/chunk_

  cat ac* | pgloader $acc_ctl
  cat refs* | pgloader $ref_ctl
  cat location* | pgloader $coord_ctl
  cat short* | pgloader $short_ctl
  cat long* | pgloader $long_ctl
  sort -u related* | pgloader $related_ctl
  cat features* | pgloader $feature_ctl

  rnac run-release
  """
}
