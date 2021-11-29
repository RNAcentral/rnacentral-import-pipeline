#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { query as locus_query } from './utils'

process merge_and_split {
  errorStrategy 'finish'

  input:
  tuple path(sequence), path(locus)

  output:
  path('by-assembly/*.json')

  """
  search-export genes select-and-split $locus $sequence by-assembly
  if [[ ! -e by-assembly/*.json ]]; then
    touch by-assembly/empty.json
  fi
  """
}

process as_xml {
  tag { "$assembly" }
  memory 10.GB
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(assembly), path('raw*.json')

  output:
  path "${xml}.gz", emit: xml
  path "count", emit: counts

  script:
  xml = "genes_${assembly}_xml.xml"
  """
  cat raw*.json > members.json
  search-export genes merge-assembly members.json merged.json
  search-export genes as-xml merged.json $xml count
  xmllint $xml --schema ${params.export.search.schema} --stream
  gzip $xml
  touch ${xml}.gz
  """
}

workflow genes {
  take:
    max_count
    sequence_json
  emit:
    xml
    counts
  main:
    Channel.fromPath('files/search-export/genes/region-info.sql') | set { locus_sql }

    locus_query(max_count, locus_sql) | set { locus_info }

    sequence_json \
    | combine(locus_info) \
    | merge_and_split \
    | flatten \
    | filter { f -> !f.isEmpty() } \
    | map { fn -> [fn.name, fn] } \
    | groupTuple \
    | as_xml

    as_xml.out.xml | set { xml }
    as_xml.out.counts | set { counts }
}
