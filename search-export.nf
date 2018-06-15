#!/usr/bin/env nextflow

process find_chunks {
  output:
  stdout raw_ranges

  """
  rnac search-export ranges ${params.search_export.chunk_size}
  """
}

ranges = raw_ranges.splitCsv()

search_query = Channel.fromPath('files/search-export/query.sql')
process export_chunk {
  publishDir params.search_export.publish, mode: 'move'
  maxForks 6

  input:
  set val(min), val(max) from ranges
  file query from search_query

  output:
  file "${xml}.gz" into search_chunks
  file("count") into search_counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  psql --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" | rnac search-export as-xml - $xml count
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
