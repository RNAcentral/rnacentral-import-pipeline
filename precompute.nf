#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { repeats } from './workflows/repeats'

process build_precompute_context {
  input:
  path('species-repeats*')

  output:
  path('context')

  """
  mkdir repeat-tree
  mv species-repeats* repeat-tree
  pushd repeat-tree
  rnac repeats build-tree .
  popd
  mkdir context
  mv repeat-tree context
  """
}

process fetch_all_urs_taxid {
  when { params.precompute.run }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query)

  output:
  path('data.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > data.csv
  """
}

process fetch_release_info {
  when { params.precompute.run }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  memory '10GB'

  input:
  path(query)

  output:
  path('data.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > raw
  sort -u raw > sorted.csv
  precompute max-release sorted.csv data.csv
  """
}

process build_urs_table {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path(load), path('xref.csv'), path('precompute.csv'), path('active.txt')

  output:
  val 'done', emit: 'flag'
  tuple path('active.txt'), path('urs.csv'), emit: to_split

  """
  precompute select xref.csv precompute.csv urs.csv
  expand-urs text active.txt urs.csv to-load.csv
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=${params.precompute.tablename} \
    -f "$load" "$PGDATABASE"
  """
}

process partial_query {
  tag { "${query.baseName}" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query)

  output:
  path("${query.baseName}.json")

  """
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=${params.precompute.tablename} \
    -f $query \
    "$PGDATABASE" > ${query.baseName}.json
  """
}

process index_data {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(data_files)

  output:
  path('index.db')

  """
  echo "${data_files.join("\n")}" > file-list
  kv index-files file-list index.db
  """
}

process build_chunks {
  input:
  tuple path('active.txt'), path('urs.txt')

  output:
  tuple path('active.txt'), path('parts/*.txt')

  script:
  def chunk_size = params.precompute.max_entries
  """
  mkdir parts
  split \
   --elide-empty-files \
   --numeric-suffixes \
   --lines=${chunk_size} \
   --additional-suffix='.txt' \
   urs.txt parts/chunk-
  """
}

process process_range {
  tag { "$urs.baseName" }
  memory params.precompute.range.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path(active), path(urs), path(context), path(db)

  output:
  path 'precompute.csv', emit: data
  path 'qa.csv', emit: qa

  """
  expand-urs text $active $urs - | kv lookup $db - raw-precompute.json
  rnac precompute from-file $context raw-precompute
  """
}

process load_data {
  beforeScript 'slack db-work loading-precompute || true'
  afterScript 'slack db-done loading-precompute || true'
  container ''

  input:
  tuple path('precompute*.csv'), path('qa*.csv'), path(pre_ctl), path(qa_ctl), path(post)

  script:
  def tablename = params.precompute.tablename
  """
  split-and-load $pre_ctl 'precompute*.csv' ${params.import_data.chunk_size} precompute
  split-and-load $qa_ctl 'qa*.csv' ${params.import_data.chunk_size} qa
  psql -v ON_ERROR_STOP=1 -f $post "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -c 'DROP TABLE IF EXISTS $tablename' "$PGDATABASE"
  """
}

// Hack to reuse the fetch_release_info process
workflow precompute_releases {
  emit: info
  main:
    Channel.fromPath('files/precompute/fetch-precompute-info.sql') \
    | fetch_release_info \
    | set { info }
}

workflow xref_releases {
  emit: info
  main:
    Channel.fromPath('files/precompute/fetch-xref-info.sql') \
    | fetch_release_info \
    | set { info }
}

workflow precompute {
  main:
    Channel.fromPath('files/precompute/queries/*.sql') | set { queries }
    Channel.fromPath('files/precompute/load-urs.sql') | set { load_sql }
    Channel.fromPath('files/all-active-urs-taxid.sql') | set { active_sql }

    build_precompute_context(repeats()) | set { context }

    precompute_releases() | set { precompute_info }
    xref_releases() | set { xref_info }
    active_sql | fetch_all_urs_taxid | set { active_urs }

    load_sql \
    | combine(xref_info) \
    | combine(precompute_info) \
    | combine(active_urs) \
    | build_urs_table

    queries \
    | combine(build_urs_table.out.flag) \
    | map { q, _ -> q } \
    | partial_query \
    | collect \
    | index_data \
    | set { indexed }

    build_urs_table.out.to_split \
    | build_chunks \
    | flatMap { active, chunks -> 
      (chunks instanceof ArrayList) ? chunks.collect { [active, it] } : [[active, chunks]]
    } \
    | filter { _, f -> !f.empty() } \
    | combine(context) \
    | combine(indexed) \
    | process_range

    process_range.out.data | collect | set { data }
    process_range.out.qa | collect | set { qa }

    data \
    | combine(qa) \
    | combine(Channel.fromPath('files/precompute/load.ctl')) \
    | combine(Channel.fromPath('files/precompute/qa.ctl')) \
    | combine(Channel.fromPath('files/precompute/post-load.sql')) \
    | load_data
}

workflow {
  precompute()
}
