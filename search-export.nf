#!/usr/bin/env nextflow

temp_publish_dir = params.search_export.publish + '-temp'
file(temp_publish_dir).mkdirs()

process find_chunks {
  output:
  file('ranges.txt') into raw_ranges

  """
  rm ${temp_publish_dir}/* || true
  rnac upi-ranges ${params.search_export.max_entries} ranges.txt
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/search-export/query.sql'))
    .into { ranges }

process export_chunk {
  publishDir temp_publish_dir, mode: 'copy'
  maxForks params.search_export.max_forks

  input:
  set val(min), val(max), file(query) from ranges

  output:
  file "${xml}.gz" into search_chunks
  file "count" into search_counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  psql --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" > search.json
  rnac search-export as-xml search.json $xml count
  xmllint $xml --schema ${params.search_export.schema} --stream
  gzip $xml
  """
}

process create_release_note {
  input:
  file('count*') from search_counts.collect()

  output:
  file 'release_note.txt' into release_note

  """
  rnac search-export release-note release_note.txt count*
  """
}

process publish {
  input:
  file 'release_note.txt' from release_note

  script:
  publish = params.search_export.publish
  """
  rm ${temp_publish_dir}/*.xml.gz
  rm ${temp_publish_dir}/release_note.txt

  cp release_note.txt ${publish}/
  mv ${temp_publish_dir}/*.xml.gz ${publish}/
  rm -r ${temp_publish_dir}
  """
}
