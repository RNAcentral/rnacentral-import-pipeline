include { query as coordinate_query} from './partial'
include { query as rfam_query} from './partial'
include { query as r2dt_query} from './partial'
include { query as prev_query} from './partial'
include { query as xref_query} from './partial'

process build {
  input:
  tuple path(coordinates), path(rfam_hits), path(r2dt_hits), path(prev), path(xref), path(chunks)

  output:
  path("chunks/*")

  """
  precompute metadata merge $coordinates $rfam_hits $r2dt_hits $prev $xref merged.json
  mkdir chunks
  precompute metadata chunk $chunks merged.json chunks
  """
}

workflow build_metadata {
  take: ranges
  emit: chunks
  main:
    Channel.fromPath('files/precompute/queries/coordinates.sql') | set { coordinates_sql }
    Channel.fromPath('files/precompute/queries/rfam-hits.sql') | set { rfam_sql }
    Channel.fromPath('files/precompute/queries/r2dt-hits.sql') | set { r2dt_sql }
    Channel.fromPath('files/precompute/queries/previous.sql') | set { prev_sql }
    Channel.fromPath('files/precompute/queries/xref.sql') | set { xref_sql }

    build(
      coordinate_query(coordinate_sql),
      rfam_query(rfam_sql),
      r2dt_query(r2dt_sql),
      prev_query(prev_sql),
      xref_query(xref_sql),
      ranges,
    ) \
    | flatten \
    | set { chunks }
}
