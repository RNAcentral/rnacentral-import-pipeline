#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { repeats } from './workflows/precompute/repeats'
include { build_precompute_accessions } from './workflows/precompute/accessions'
include { build_urs_table } from './workflows/precompute/build_urs_table'

include { query as coordinate_query} from './workflows/precompute/utils'
include { query as rfam_query} from './workflows/precompute/utils'
include { query as r2dt_query} from './workflows/precompute/utils'
include { query as prev_query} from './workflows/precompute/utils'
include { query as xref_query} from './workflows/precompute/utils'

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

process build_metadata {
  input:
  tuple path(coordinates), path(rfam_hits), path(r2dt_hits), path(prev), path(xref)

  output:
  path("metadata.json")

  """
  precompute metadata merge $coordinates $rfam_hits $r2dt_hits $prev $xref metadata.json
  """
}

process build_ranges {
  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.precompute.max_entries
  def tablename = params.precompute.tablename
  """
  rnac upi-ranges --table-name $tablename $chunk_size ranges.csv
  """
}

process query_accession_range {
  tag { "$min-$max" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  maxForks params.precompute.maxForks

  input:
  tuple val(min), val(max), path(query)

  output:
  tuple val(min), val(max), path('accessions.json')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v min=$min \
    -v max=$max \
    -f $query \
    "$PGDATABASE" > accessions.json
  """
}

process process_range {
  tag { "$min-$max" }
  memory params.precompute.range.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(accessions), path(metadata), path(context)

  output:
  path 'precompute.csv', emit: data
  path 'qa.csv', emit: qa

  """
  precompute normalize $accessions $metadata merged.json
  rnac precompute from-file $context merged.json
  """
}

process load_data {
  beforeScript 'slack db-work loading-precompute || true'
  afterScript 'slack db-done loading-precompute || true'
  container ''

  input:
  path('precompute*.csv') 
  path('qa*.csv')
  path(pre_ctl)
  path(qa_ctl)
  path(post)

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
  main:
    Channel.fromPath('files/precompute/get-accessions/query.sql') | set { accession_query }
    Channel.fromPath('files/precompute/load.ctl') | set { data_ctl }
    Channel.fromPath('files/precompute/qa.ctl') | set { qa_ctl }
    Channel.fromPath('files/precompute/post-load.sql') | set { post_load }

    // Various metadata components for precompute
    Channel.fromPath('files/precompute/queries/coordinates.sql') | set { coordinate_sql }
    Channel.fromPath('files/precompute/queries/rfam-hits.sql') | set { rfam_sql }
    Channel.fromPath('files/precompute/queries/r2dt-hits.sql') | set { r2dt_sql }
    Channel.fromPath('files/precompute/queries/previous.sql') | set { prev_sql }
    Channel.fromPath('files/precompute/queries/xref.sql') | set { xref_sql }

    repeats | build_precompute_context | set { context }
    build_urs_table | set { built_table }
    built_table | build_precompute_accessions | set { accessions_ready }

    build_metadata(
      coordinate_query(built_table, coordinate_sql),
      rfam_query(built_table, rfam_sql),
      r2dt_query(built_table, r2dt_sql),
      prev_query(built_table, prev_sql),
      xref_query(built_table, xref_sql),
    ) \
    | set { metadata }

    built_table \
    | build_ranges \
    | splitCsv \
    | combine(accessions_ready) \
    | combine(accession_query) \
    | map { _tablename, min, max, _flag, sql -> [min, max, sql] } \
    | query_accession_range \
    | combine(metadata) \
    | combine(context) \
    | process_range

    process_range.out.data | collect | set { data }
    process_range.out.qa | collect | set { qa }

    load_data(data, qa, data_ctl, qa_ctl, post_load)
}

workflow {
  precompute()
}
