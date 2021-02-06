#!/usr/bin/env nextflow

include { query as accession_query } from './workflows/search-export/utils'
include { query as base_query } from './workflows/search-export/utils'
include { query as crs_query } from './workflows/search-export/utils'
include { query as feeback_query } from './workflows/search-export/utils'
include { query as go_query } from './workflows/search-export/utils'
include { query as prot_query } from './workflows/search-export/utils'
include { query as rnas_query } from './workflows/search-export/utils'
include { query as precompute_query } from './workflows/search-export/utils'
include { query as qa_query } from './workflows/search-export/utils'
include { query as r2dt_query } from './workflows/search-export/utils'
include { query as ref_query } from './workflows/search-export/utils'
include { query as rfam_query } from './workflows/search-export/utils'

process fetch_so_tree {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query)

  output:
  path('so-term-tree.json')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  rnac search-export so-term-tree raw.json so-term-tree.json
  """
}

process build_json {
  input:
  path(base)
  path(crs)
  path(feeback)
  path(go)
  path(prot)
  path(rnas)
  path(precompute)
  path(qa)
  path(r2dt)
  path(rfam)
  path(so_tree)

  output:
  path("merged.json")

  """
  search-export merge $base $crs $feeback $go $prot $rnas $precompute $qa $r2dt $rfam $so_tree merged.json
  """
}

process build_ranges {
  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.search_export.max_entries
  """
  rnac upi-ranges --table-name search_export_accessions $chunk_size ranges.csv
  """
}
 
process fetch_accession {
  tag { "$min-$max" }

  input:
  tuple val(min), val(max), val(sql)

  output:
  path("raw.json")

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v min-$min \
    -v max=$max \
    -f "$sql" \
    "$PGDATABASE" > raw.json
  """
}

process export_chunk {
  tag { "$index" }
  memory params.search_export.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(accessions), path(metadata)

  output:
  path "${xml}.gz", emit: xml
  path "count", emit: counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  search-export normalize $accessions $metadata raw.json
  rnac search-export as-xml raw.json $xml count
  xmllint $xml --schema ${params.search_export.schema} --stream
  gzip $xml
  """
}

process create_release_note {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('count*')

  output:
  path('release_note.txt')

  script:
  """
  rnac search-export release-note ${params.release} release_note.txt count*
  """
}

// At this point we should be able to safely move data into the final location.
// This deletes the old data and then moves the new data in place.
process atomic_publish {
  container ''

  input:
  path('release_note.txt')
  path(xml)

  script:
  def publish = params.search_export.publish
  if (params.search_export.publish.host)
    """
    ssh "$publish.host" 'mkdir -p $publish.path' || true
    ssh "$publish.host" 'rm -r $publish.path/*' || true
    scp ${xml} 'release_note.txt' $publish.host:$publish.path
    """
  else
    """
    rm $publish.path/*
    cp ${xml} 'release_note.txt' $publish.path
    """
}

workflow search_export {
  Channel.fromPath('files/search-export/parts/base.sql') | set { base_sql }
  Channel.fromPath('files/search-export/parts/crs.sql') | set { crs_sql }
  Channel.fromPath('files/search-export/parts/feedback.sql') | set { feeback_sql }
  Channel.fromPath('files/search-export/parts/go-annotations.sql') | set { go_sql }
  Channel.fromPath('files/search-export/parts/interacting-proteins.sql') | set { prot_sql }
  Channel.fromPath('files/search-export/parts/interacting-rnas.sql') | set { rnas_sql }
  Channel.fromPath('files/search-export/parts/precompute.sql') | set { precompute_sql }
  Channel.fromPath('files/search-export/parts/qa-status.sql') | set { qa_sql }
  Channel.fromPath('files/search-export/parts/r2dt.sql') | set { r2dt_sql }
  Channel.fromPath('files/search-export/parts/rfam-hits.sql') | set { rfam_sql }

  Channel.fromPath('files/search-export/get-accessions.sql') | set { accession_sql }

  build_search_table | set { search_ready }
  search_ready | build_search_accessions | set { accessions_ready }

  build_json(
    base_query(search_ready, base_sql),
    crs_query(search_ready, crs_sql),
    feeback_query(search_ready, feeback_sql),
    go_query(search_ready, go_sql),
    prot_query(search_ready, prot_sql),
    rnas_query(search_ready, rnas_sql),
    precompute_query(search_ready, precompute_sql),
    qa_query(search_ready, qa_sql),
    r2dt_query(search_ready, r2dt_sql),
    rfam_query(search_ready, rfam_sql),
    fetch_so_tree(search_ready, so_sql),
  )\
  | set { metadata }

  search_ready \
  | build_ranges \
  | map { _tablename, min, max -> [min, max ] } \
  | combine(accessions_sql) \
  | fetch_accession \
  | combine(metadata) \
  | export_chunk

  export_chunk.out.counts | collect | create_release_note | set { note }
  export_chunk.out.xml | collect | set { xml }

  atomic_publish(note, xml)
}

workflow {
  search_export()
}
