#!/usr/bin/env nextflow

process find_chunks {
  output:
  file('ranges.txt') into raw_ranges

  """
  rnac upi-ranges ${params.search_export.max_entries} ranges.txt
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/search-export/query.sql'))
  .into { ranges }

process export_chunk {
  publishDir params.search_export.publish, mode: 'copy'
  maxForks 6

  input:
  set val(min), val(max), file(query) from ranges

  output:
  file "${xml}.gz" into search_chunks
  file("count") into search_counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  psql --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" |\
    rnac search-export as-xml - $xml count
  xmllint $xml --schema ${params.search_export.schema} --stream
  gzip $xml
  """
}

process create_release_note {
  publishDir params.search_export.publish, mode: 'move'

  input:
  file('count*') from search_counts.collect()

  output:
  file 'release_note.txt' into release_note

  script:
  """
  rnac search-export release-note release_note.txt count*
  """
}
