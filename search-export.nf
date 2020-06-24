#!/usr/bin/env nextflow

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

Channel.fromPath('files/search-export/metadata/*.sql')
  .set { metadata_queries }

process fetch_metdata {
  maxForks 2

  input:
  file(query) from metadata_queries

  output:
  file("${query.baseName}.json") into standard_metadata

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > "${query.baseName}.json"
  """
}

process fetch_so_tree {
  input:
  file(query) from Channel.fromPath('files/search-export/so-rna-types.sql')

  output:
  file('so-term-tree.json') into so_term_tree_metadata

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  rnac search-export so-term-tree raw.json so-term-tree.json
  """
}

standard_metadata
  .collect()
  .set { unmerged_metdata }

process merge_metadata {
  input:
  file(metadata) from unmerged_metdata
  file(so_data) from so_term_tree_metadata

  output:
  file('merged.db') into metadata

  """
  cat $metadata > metadata.json
  rnac search-export merge-metadata metadata.json $so_data merged.db
  """
}

process export_search_json {
  // beforeScript 'slack db-work search-export || true'
  // afterScript 'slack db-done search-export || true'
  tag { "$min-$max" }
  maxForks params.search_export.max_forks
  echo true

  input:
  set val(tablename), val(min), val(max), file(query) from ranges

  output:
  set val(min), val(max), file('search.json') into raw_json

  script:
  """
  psql -v ON_ERROR_STOP=1 --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" > search.json
  """
}

raw_json
  .filter { min, max, fn -> !fn.empty() }
  .combine(metadata)
  .set { search_json }

process export_chunk {
  tag { "$min-$max" }
  memory params.search_export.memory

  input:
  set val(min), val(max), file(json), file(metadata) from search_json

  output:
  file("${xml}.gz") into search_chunks
  file "count" into search_counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  rnac search-export as-xml ${json} ${metadata} ${xml} count
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
  rnac search-export release-note ${params.release} release_note.txt count*
  """
}

// At this point we should be able to safely move data into the final location.
// This deletes the old data and then moves the new data in place.
process atomic_publish {
  input:
  file('release_note.txt') from release_note
  file(xml) from search_chunks.collect()

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
