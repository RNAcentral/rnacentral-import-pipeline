#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { query as accession_query } from './utils'
include { query as base_query } from './utils'
include { query as crs_query } from './utils'
include { query as feeback_query } from './utils'
include { query as go_query } from './utils'
include { query as prot_query } from './utils'
include { query as rnas_query } from './utils'
include { query as precompute_query } from './utils'
include { query as qa_query } from './utils'
include { query as r2dt_query } from './utils'
include { query as ref_query } from './utils'
include { query as rfam_query } from './utils'
include { query as orf_query } from './utils'
include { query as locus_query } from './utils'
include { build_search_accessions } from './build-accession-table'

process setup {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(sql)
  path(counts)

  output:
  path('counts.txt')

  """
  psql -v ON_ERROR_STOP=1 -f "$sql" "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -f "$counts" "$PGDATABASE" > counts.txt
  """
}

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

process build_metadata {
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
  path(orf)
  path(text)
  path(so_tree)

  output:
  path("merged.json")

  """
  search-export sequences merge $base $crs $feeback $go $prot $rnas $precompute $qa $r2dt $rfam $orf $text $so_tree merged.json
  """
}

process build_ranges {
  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.export.search.max_entries
  """
  rnac upi-ranges --table-name search_export_urs $chunk_size ranges.csv
  """
}

process fetch_accession {
  tag { "$min-$max" }
  maxForks 3
  time '10m'
  errorStrategy 'retry'
  maxRetries 5

  input:
  tuple val(min), val(max), path(sql), val(_flag)

  output:
  tuple val(min), val(max), path("raw.json")

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v min=$min \
    -v max=$max \
    -f "$sql" \
    "$PGDATABASE" > raw.json
  """
}

process text_mining_query {
  input:
  val(max_count)
  path(script)

  output:
  path("publication-count.json")

  """
  curl "$params.export.search.text_mining" > counts.csv
  psql -v ON_ERROR_STOP=1 -c "\\copy search_export_publication_counts from 'counts.csv'" "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -f "$script" "$PGDATABASE" > raw.json
  search-export group publication-count raw.json ${max_count} publication-count.json
  """
}

process as_xml {
  tag { "$min-$max" }
  memory params.export.search.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(raw), path(metadata)

  output:
  path "${xml}.gz", emit: xml
  path "count", emit: counts
  path "data.json", emit: sequences

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  search-export sequences normalize $raw $metadata data.json
  rnac search-export as-xml data.json $xml count
  xmllint $xml --schema ${params.export.search.schema} --stream
  gzip $xml
  """
}

workflow sequences {
  emit:
    xml
    counts
    search_count
    sequence_json
  main:
    Channel.fromPath('files/search-export/setup.sql') | set { setup_sql }

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
    Channel.fromPath('files/search-export/parts/orfs.sql') | set { orf_sql }
    Channel.fromPath('files/search-export/parts/text-mining.sql') | set { text_sql }
    Channel.fromPath('files/search-export/so-rna-types.sql') | set { so_sql }

    Channel.fromPath('files/search-export/parts/accessions.sql') | set { accessions_sql }

    Channel.fromPath('files/search-export/get-counts.sql') | set { counts_sql }

    setup(setup_sql, counts_sql)
    | splitCsv \
    | first \
    | map { row -> row[0].toInteger() + 1 } \
    | set { search_count }

    search_count | build_search_accessions | set { accessions_ready }

    build_metadata(
      base_query(search_count, base_sql),
      crs_query(search_count, crs_sql),
      feeback_query(search_count, feeback_sql),
      go_query(search_count, go_sql),
      prot_query(search_count, prot_sql),
      rnas_query(search_count, rnas_sql),
      precompute_query(search_count, precompute_sql),
      qa_query(search_count, qa_sql),
      r2dt_query(search_count, r2dt_sql),
      rfam_query(search_count, rfam_sql),
      orf_query(search_count, orf_sql),
      text_mining_query(search_count, text_sql),
      fetch_so_tree(so_sql),
    )\
    | set { metadata }

    search_count \
    | build_ranges \
    | splitCsv \
    | map { _tablename, min, max -> [min.toInteger(), max.toInteger()] }
    | map { min, max -> [min, max - 1] } \
    | combine(accessions_sql) \
    | combine(accessions_ready) \
    | fetch_accession \
    | combine(metadata) \
    | as_xml

    as_xml.out.sequences | set { sequence_json }
    as_xml.out.counts | set { counts }
    as_xml.out.xml | set { xml }
}
