#!/usr/bin/env nextflow

queries = [
  all: 'files/precompute/methods/all.sql',
  using_release: 'files/precompute/methods/using-release.sql'
]

if (params.precompute.methods.all) {
  Channel
    .from([[file(queries.all), []]])
    .set { precompute_upi_queries }
} else {
  Channel
    .from([[file(queries.using_release), []]])
    .set { precompute_upi_queries }
}

process query_upis {
  input:
  set file(sql), val(variables) from precompute_upi_queries

  output:
  file('ranges.txt') into raw_ranges

  script:
  def tablename = params.precompute.tablename
  def variable_option = variables.inject([]) { acc, entry ->
    acc << """-v ${entry.key}="${entry.value}" """
  }
  """
  psql -v ON_ERROR_STOP=1 -f "$sql" ${variable_option.join(' ')} "$PGDATABASE"
  rnac upi-ranges --table-name $tablename ${params.precompute.max_entries} ranges.txt
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/precompute/query.sql'))
  .set { ranges }

process precompute_range_query {
  beforeScript 'slack db-work precompute-range || true'
  afterScript 'slack db-done precompute-range || true'
  maxForks params.precompute.maxForks

  input:
  set val(tablename), val(min), val(max), file(query) from ranges

  output:
  file 'raw-precompute.json' into precompute_raw

  """
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=$tablename \
    --variable min=$min \
    --variable max=$max \
    -f "$query" \
    '$PGDATABASE' > raw-precompute.json
  """
}

process precompute_range {
  memory '4 GB'

  input:
  file(raw) from precompute_raw

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


  script:
  def tablename = params.precompute.tablename
  """
  cp $pre_ctl _pre.ctl
  cp $qa_ctl _qa.ctl
  pgloader _pre.ctl
  pgloader _qa.ctl
  psql -v tablename=$tablename -v ON_ERROR_STOP=1 -f $post "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -c 'DROP TABLE IF EXISTS $tablename' "$PGDATABASE"
  """
}

process mods_for_feedback {
  input:
  file(query) from Channel.fromPath('files/precompute/find-mod-info.sql')

  output:
  file('info') into raw_mods

  script:
  names = []
  for (name in params.precompute.feedback.databases) {
    names << "'${name.toUpperCase()}'"
  }
  names = '(' + names.join(', ') + ')'
  """
  psql -v ON_ERROR_STOP=1 -v "names=${names}" -f "$query" "$PGDATABASE" > info
  cat info
  """
}

raw_mods
  .splitCsv()
  .combine(Channel.fromPath('files/ftp-export/genome_coordinates/query.sql'))
  .set { mods }

process generate_feedback_report {
  memory '4 GB'

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
  cp $ctl local_$ctl
  cat feedback*.tsv > merged.tsv
  pgloader --on-error-stop local_$ctl
  """
}
