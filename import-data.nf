#!/usr/bin/env nextflow

// This channel will store all the outputs for all data to import. We then split
// the channel into several different channels at the end. Doing this prevents
// us from having many large duplicated blocks. Also adding a new file type to
// import is much easier this way, as we only have to modify in a few places
// instead of every possible process that could produce it.
raw_output = Channel.empty()
dataless_imports = Channel.empty()
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
  } else if (database.containsKey('remote') && database.remote) {
    db_channel = Channel.value([name, database])
    to_fetch.mix(db_channel).set { to_fetch }
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

process fetch_data {
  input:
  set val(name), val(database) from to_fetch

  output:
  set val(name), file("${pattern}") into fetched

  script:
  pattern = database.pattern
  find_cmd = []
  for (pattern in database.get('excluded_patterns', [])) {
    find_cmd << "-name '${pattern}'"
  }

  if (find_cmd.size) {
    find_cmd = "find . ${find_cmd.join(' ')} | xargs rm"
  } else {
    find_cmd = ''
  }

  """
  fetch ${database.remote} $pattern
  ${find_cmd}
  """
}

metadataless = Channel.create()
with_metadata = Channel.create()
fetched
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
  if (filename.toString().endsWith(".gz")) {
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
  when:
  params.import_data.databases['rfam'] || params.import_data.databases['ensembl']

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
  when:
  params.import_data.databases['ena']

  input:
  file tpa_file from tpa_url_file

  output:
  set val('ena'), file("tpa.tsv") into ena_metadata

  """
  cat $tpa_file | xargs wget -O - >> tpa.tsv
  """

}

process fetch_rgd_metadata {
  when:
  params.import_data.databases['rgd']

  output:
  set val('rgd'), file('genes.txt') into rgd_metadata

  """
  wget ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt > genes.txt
  """
}

Channel.empty()
  .mix(ena_metadata)
  .mix(ensembl_metadata)
  .mix(rfam_metadata)
  .mix(rgd_metadata)
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
