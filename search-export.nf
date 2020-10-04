#!/usr/bin/env nextflow

process find_chunks {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query) from Channel.fromPath('files/find-active-xrefs-urs-taxids.sql')

  output:
  path('parts/*') into ids

  script:
  """
  mkdir parts
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" | sort -u > upi-taxids
  split -n ${params.search_export.max_entries} upi-taxids parts/
  """
}

process fetch_parts {
  maxForks 2
  container ''

  input:
  path(query) from query_parts

  output:
  path("${query.baseName}.json") into Channel.fromPath('files/search-export/parts/*.sql')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > "${query.baseName}.json"
  """
}

process fetch_so_tree {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query) from Channel.fromPath('files/search-export/so-rna-types.sql')

  output:
  path('so-term-tree.json') into so_term_tree_metadata

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > raw.json
  rnac search-export so-tree raw.json so-term-tree.json
  """
}

process index_parts {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('*.json') from query_parts.collect()

  output:
  path('merged.db') into index

  """
  find . -name '*.json' | kv index-files - merged.db
  """
}

ids
  .combine(so_term_tree_metadata)
  .combine(index)
  .set { search_data }

process export_chunk {
  tag { "$min-$max" }
  memory params.search_export.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path(id_file), path(so_tree), path(database) from search_data

  output:
  path("${xml}.gz") into search_chunks
  path("count") into search_counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  kv fetch $id_file $database - | search-utils normalize - $so_tree - | rnac search-export as-xml - ${xml} count
  xmllint ${xml} --schema ${params.search_export.schema} --stream
  gzip ${xml}
  """
}

process create_release_note {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('count*') from search_counts.collect()

  output:
  path('release_note.txt') into release_note

  """
  rnac search-export release-note ${params.release} release_note.txt count*
  """
}

// At this point we should be able to safely move data into the final location.
// This deletes the old data and then moves the new data in place.
process atomic_publish {
  container ''

  input:
  path('release_note.txt') from release_note
  path(xml) from search_chunks.collect()

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
