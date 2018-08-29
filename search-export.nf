#!/usr/bin/env nextflow

// This is used to provide a temporary directory publish to. We publish here
// and then at the we publish everything at once. This prevents writing part of
// the data in the case one job fails while already succeeded.
tmp = "mktemp -d -p ${workDir}".execute().text.trim()

process find_chunks {
  output:
  file('ranges.txt') into raw_ranges

  script:
  """
  rnac upi-ranges ${params.search_export.max_entries} ranges.txt
  """
}

raw_ranges
  .splitCsv()
  .combine(Channel.fromPath('files/search-export/query.sql'))
    .set { ranges }

process export_search_json {
  // beforeScript 'slack db-work search-export || true'
  // afterScript 'slack db-done search-export || true'
  maxForks params.search_export.max_forks

  input:
  set val(min), val(max), file(query) from ranges

  output:
  set val(min), val(max), file('search.json') into search_json

  script:
  """
  psql --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" > search.json
  """
}

process export_chunk {
  publishDir "${tmp}/", mode: 'copy'

  input:
  set val(min), val(max), file(json) from search_json

  output:
  file("${xml}.gz") into search_chunks
  file "count" into search_counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  rnac search-export as-xml ${json} ${xml} count
  xmllint ${xml} --schema ${params.search_export.schema} --stream
  gzip ${xml}
  """
}

process create_release_note {
  input:
  file('count*') from search_counts.collect()

  output:
  file 'release_note.txt' into release_note

  script:
  """
  rnac search-export release-note release_note.txt count*
  """
}

// At this point we should be able to safely move data into the final location.
// This deletes the old data and then moves the new data in place.
process atomic_publish {
  input:
  file(release) from release_note

  script:
  dir = params.search_export.publish
  """
  [ -d $dir ] || mkdir -p $dir
  rm $dir/*.xml.gz || true

  mv ${tmp}/*.xml.gz $dir/
  cp $release $dir/release_note.txt
  """
}
