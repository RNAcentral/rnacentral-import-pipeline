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
  // import. For example if custom was ./custom we would expect it to contain:
  // ./custom/ac_info.csv, ./custom/genomic_locations.csv, etc. The files may be
  // nested so: ./custom/mirbase/ac_info.csv would also work.
  if (name == "custom") {
    raw_output.mix(Channel.fromPath("${base}/**.csv")).flatten().set { raw_output }
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

  script:
  // wget doesn't play well with ** so we use lftp instead. Also, wget doesn't
  // like following symlinks, but lftp does. I should find a way to indicate
  // this to the pipeline.
  if (remote.startsWith('ftp:') && remote.contains('**')) {
    """
    lftp -c 'mget ${remote}'
    find . -name '*.tar.gz' | xargs tar xvf
    """
  } else {
    """
    wget "${remote}"
    find . -name '*.tar.gz' | xargs tar xvf
    """
  }
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
  memory params.databases[name].get('memory', '2 GB')

  input:
  set val(name), file(filename) from metadataless

  output:
  file "*.csv" into raw_metadataless_output mode flatten

  script:
  if (filename.endsWith(".gz")) {
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
  memory params.databases[name].get('memory', '2 GB')

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

raw_output
  .map { f ->
    filename = f.getName()
    [filename.take(filename.lastIndexOf('.')), f]
  }
  .groupTuple()
  .set { to_merge }

process merge_csvs {
  input:
  val(name), file('raw*.csv') from to_merge

  output:
  file("merged/${name}*.csv") into merged mode flatten

  """
  find . -name 'raw*.csv' |\
  xargs cat |\
  split --additional-suffix=.csv -dC ${params.import_data.chunk_size} - merged/${name}
  """
}

accessions = Channel.create()
locations = Channel.create()
long_sequences = Channel.create()
short_sequences = Channel.create()
refs = Channel.create()
secondary = Channel.create()
related = Channel.create()
features = Channel.create()

merged
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
  file 'ac_info*.csv' from accessions.collect()
  file 'locations*.csv' from locations.collect()
  file 'long*.csv' from long_sequences.collect()
  file 'short*.csv' from short_sequences.collect()
  file 'refs*.csv' from refs.collect()
  file 'related*.csv' from related.collect()
  file 'features*.csv' from features.collect()

  file acc_ctl from Channel.fromPath('files/import-data/accessions.ctl')
  file locations_ctl from Channel.fromPath('files/import-data/locations.ctl')
  file long_ctl from Channel.fromPath('files/import-data/long-sequences.ctl')
  file short_ctl from Channel.fromPath('files/import-data/short-sequences.ctl')
  file ref_ctl from Channel.fromPath('files/import-data/references.ctl')
  file related_ctl from Channel.fromPath('files/import-data/related-sequences.ctl')
  file features_ctl from Channel.fromPath('files/import-data/features.ctl')

  """
  find . -name '.ctl' | xargs -I {} cp {} local_{}

  pgloader local_${acc_ctl}
  pgloader local_${locations_ctl}
  pgloader local_${long.ctl}
  pgloader local_${short.ctl}
  pgloader local_${ref_ctl}

  rnac run-release

  pgloader local_${features_ctl}
  sort -u related* | pgloader local_$related_ctl
  """
}
