#!/usr/bin/env nextflow

process find_chunks {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  output:
  path('ranges.txt')

  """
  rnac upi-ranges ${params.search_export.max_entries} ranges.txt
  """
}

process fetch_metdata {
  maxForks 2
  container ''

  input:
  path(query)

  output:
  path("${query.baseName}.json")

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > "${query.baseName}.json"
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

process merge_metadata {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(metadata)
  path(so_data)

  output:
  path('merged.db')

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
  container ''

  input:
  tuple val(tablename), val(min), val(max), path(query)

  output:
  tuple val(min), val(max), path('search.json')

  script:
  """
  psql -v ON_ERROR_STOP=1 --variable min=$min --variable max=$max -f "$query" "$PGDATABASE" > search.json
  """
}

process export_chunk {
  tag { "$min-$max" }
  memory params.search_export.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(json), path(metadata)

  output:
  path "${xml}.gz", emit: xml
  path "count", emit: counts

  script:
  xml = "xml4dbdumps__${min}__${max}.xml"
  """
  rnac search-export as-xml ${json} ${metadata} ${xml} count
  xmllint ${xml} --schema ${params.search_export.schema} --stream
  gzip ${xml}
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
  Channel.fromPath('files/search-export/query.sql').set { query_sql}
  Channel.fromPath('files/search-export/metadata/*.sql').set { metadata_queries }
  Channel.fromPath('files/search-export/so-rna-types.sql').set { so_query }

  so_query | fetch_so_tree | set { so_tree }
  metadata_queries | fetch_metadata | collect | set { raw_metadata }
  merge_metadata(raw_metadata, so_tree) | set { metadata }

  find_chunks \
  | splitCsv \
  | combine(query_sql) \
  | export_search_json \
  | filter { min, max, fn -> !fn.empty() } \
  | combine(metadata) \
  | export_chunk

  export_chunk.out.counts | collect | create_release_note | set { note }
  export_chunk.out.xml | collect | set { xml }

  atomic_publish(note, xml)
}

workflow {
  search_export()
}
