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

process find_ranges {
  when { params.precompute.run }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(sql)

  output:
  path('ranges.txt')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=$params.precompute.tablename \
    -f "$sql" "$PGDATABASE"
  rnac upi-ranges --table-name $params.precompute.tablename ${params.precompute.max_entries} ranges.txt
  """
}

process query_range {
  tag { "$min-$max" }
  beforeScript 'slack db-work precompute-range || true'
  afterScript 'slack db-done precompute-range || true'
  maxForks params.precompute.maxForks
  container ''

  input:
  tuple val(tablename), val(min), val(max), path(query)

  output:
  tuple val(min), val(max), path('raw-precompute.json')

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

process process_range {
  tag { "$min-$max" }
  memory params.precompute.range.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(raw), path(context)

  output:
  path('precompute.csv'), emit: data
  path('qa.csv'), emit: qa

  """
  rnac precompute from-file $context $raw
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

workflow precompute {
  take: method

  main:
    build_precompute_context(repeats()) | set { context }

    method \
    | find_ranges \
    | splitCsv \
    | combine(Channel.fromPath('files/precompute/query.sql')) \
    | query_range \
    | combine(context) \
    | process_range \

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
  precompute(Channel.fromPath("files/precompute/methods/${params.precompute.method.replace('_', '-')}.sql"))
}
