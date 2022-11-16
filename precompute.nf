#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* include { repeats } from './workflows/precompute/repeats' */
include { build_precompute_accessions } from './workflows/precompute/accessions'
include { build_urs_table } from './workflows/precompute/build_urs_table'

include { query as coordinate_query} from './workflows/precompute/utils'
include { query as rfam_query} from './workflows/precompute/utils'
include { query as r2dt_query} from './workflows/precompute/utils'
include { query as prev_query} from './workflows/precompute/utils'
include { query as basic_query} from './workflows/precompute/utils'
include { query as orf_query} from './workflows/precompute/utils'

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

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
  path(basic)
  path(coordinates)
  path(rfam_hits)
  path(r2dt_hits)
  path(prev)
  path(orf)

  output:
  path("metadata.json")

  """
  precompute metadata merge $basic $coordinates $rfam_hits $r2dt_hits $prev $orf metadata.json
  """
}

process build_ranges {
  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.precompute.max_entries
  """
  rnac upi-ranges --table-name precompute_urs $chunk_size ranges.csv
  """
}

process find_upi_taxid_ranges {
  executor 'local'

  input:
  path(ranges)

  output:
  path('urs_taxid.csv')

  """
  rnac precompute upi-taxid-ranges $ranges urs_taxid.csv
  """
}

process query_accession_range {
  tag { "$min-$max" }
  maxForks params.precompute.maxForks
  memory '4GB'

  input:
  tuple val(min), val(max), path(query), val(upi_start), val(upi_stop)

  output:
  tuple val(min), val(max), path('accessions.json')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v min=$min \
    -v max=$max \
    -f $query \
    "$PGDATABASE" > raw.json
  precompute group-accessions raw.json $upi_start $upi_stop accessions.json
  """
}

process process_range {
  tag { "$min-$max" }
  memory params.precompute.range.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(accessions), path(metadata)

  output:
  path 'precompute.csv', emit: data
  path 'qa.csv', emit: qa

  """
  mkdir context
  precompute normalize $accessions $metadata merged.json
  rnac precompute from-file context merged.json
  """
}

process load_data {

  input:
  path('precompute*.csv')
  path('qa*.csv')
  path(pre_ctl)
  path(qa_ctl)
  path(post)

  script:
  """
  split-and-load $pre_ctl 'precompute*.csv' ${params.import_data.chunk_size} precompute
  split-and-load $qa_ctl 'qa*.csv' ${params.import_data.chunk_size} qa
  psql -v ON_ERROR_STOP=1 -f $post "$PGDATABASE"
  """
}

workflow precompute {
  take: _flag
  main:

    Channel.of("Starting precompute pipeline") | slack_message

    Channel.fromPath('files/precompute/get-accessions/query.sql') | set { accession_query }
    Channel.fromPath('files/precompute/load.ctl') | set { data_ctl }
    Channel.fromPath('files/precompute/qa.ctl') | set { qa_ctl }
    Channel.fromPath('files/precompute/post-load.sql') | set { post_load }

    Channel.fromPath('files/precompute/queries/basic.sql') | set { basic_sql }
    Channel.fromPath('files/precompute/queries/coordinates.sql') | set { coordinate_sql }
    Channel.fromPath('files/precompute/queries/rfam-hits.sql') | set { rfam_sql }
    Channel.fromPath('files/precompute/queries/r2dt-hits.sql') | set { r2dt_sql }
    Channel.fromPath('files/precompute/queries/previous.sql') | set { prev_sql }
    Channel.fromPath('files/precompute/queries/orfs.sql') | set { orf_sql }

    // repeats | build_precompute_context | set { context }
    Channel.of(params.precompute.method) | build_urs_table | set { urs_counts }
    urs_counts | build_precompute_accessions | set { accessions_ready }

    build_metadata(
      basic_query(urs_counts, basic_sql),
      coordinate_query(urs_counts, coordinate_sql),
      rfam_query(urs_counts, rfam_sql),
      r2dt_query(urs_counts, r2dt_sql),
      prev_query(urs_counts, prev_sql),
      orf_query(urs_counts, orf_sql),
    ) \
    | set { metadata }

    urs_counts \
    | build_ranges \
    | set { ranges }

    ranges \
    | find_upi_taxid_ranges \
    | splitCsv \
    | map { upi_min, ut_min, ut_max -> [upi_min, ut_min, ut_max.toInteger() + 1] } \
    | set { upi_taxid_ranges }

    ranges \
    | splitCsv \
    | combine(accessions_ready) \
    | combine(accession_query) \
    | map { _tablename, min, max, _flag, sql -> [min, max, sql] } \
    | join(upi_taxid_ranges) \
    | query_accession_range \
    | combine(metadata) \
    | process_range

    process_range.out.data | collect | set { data }
    process_range.out.qa | collect | set { qa }

    load_data(data, qa, data_ctl, qa_ctl, post_load)
}

workflow {
  precompute(Channel.of(true))
}

workflow.onComplete {

  slack_closure("Precompute workflow completed. Data import complete")
}

workflow.onError {

  slack_closure("Precompute workflow encountered an error and crashed")
}
