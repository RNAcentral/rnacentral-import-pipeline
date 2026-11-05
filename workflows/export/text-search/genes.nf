#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { query as locus_query } from './utils'
include { fetch_schema } from './utils'

process merge_and_split {
  // SLURM kills carry no exit status, so retry on the attempt count and grow memory.
  memory { 2.GB * task.attempt }
  errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
  maxRetries 3

  input:
  tuple path(sequence), path(locus)

  output:
  path('by-assembly/*.json')

  script:
  """
  search-export genes select-and-split $locus $sequence by-assembly
  mkdir -p by-assembly
  compgen -G 'by-assembly/*.json' > /dev/null || touch by-assembly/empty.json
  """
}

process as_xml {
  tag { "$assembly" }
  // SLURM kills carry no exit status, so retry on the attempt count and grow memory.
  memory { params.export.search.memory * task.attempt }
  errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
  maxRetries 3

  input:
  tuple val(assembly), path('raw*.json'), path('so_tree.json'), path('schema.xsd')

  output:
  path "${xml}.gz", emit: xml
  path "count", emit: counts

  script:
  xml = "genes_${assembly}_xml.xml"
  """
  cat raw*.json > members.json
  search-export genes merge-assembly so_tree.json members.json merged.json
  search-export genes as-xml merged.json $xml count
  xmllint $xml --schema schema.xsd --stream
  gzip $xml
  """
}

workflow genes {
  take:
    max_count
    sequence_json
    so_tree
  main:
    channel.fromPath('files/search-export/genes/region-info.sql') | set { locus_sql }

    locus_query(max_count, locus_sql) | set { locus_info }
    fetch_schema()

    sequence_json \
    | combine(locus_info) \
    | merge_and_split \
    | flatten \
    | filter { f -> !f.isEmpty() } \
    | map { fn -> [fn.name, fn] } \
    | groupTuple \
    | combine(so_tree) \
    | combine(fetch_schema.out) \
    | as_xml

    as_xml.out.xml | set { xml }
    as_xml.out.counts | set { counts }
  emit:
    xml
    counts
}
