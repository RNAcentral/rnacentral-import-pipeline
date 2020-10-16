#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_taxids {
  input:
  path(query)

  output:
  path('parts/*')

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
  path(query)

  output:
  path("${query.baseName}.json") into Channel.fromPath('files/search-export/parts/*.sql')

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
  rnac search-export so-tree raw.json so-term-tree.json
  """
}

process index_parts {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('*.json')

  output:
  path('merged.db')

  """
  find . -name '*.json' | kv index-files - merged.db
  """
}

process export_chunk {
  tag { "$min-$max" }
  memory params.search_export.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path(id_file), path(database), path(so_tree)

  output:
  path("${xml}.gz") emit: xml
  path("count") emit: counts

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
  path('count*')

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
  Channel.fromPath('files/search-export/so-rna-types.sql') \
  | fetch_so_tree \
  | set { so_tree }

  Channel.fromPath('files/search-export/parts/*.sql') \
  | fetch_parts \
  | collect \
  | index_parts \
  | combine(so_tree) \
  | set { data }

  Channel.fromPath('files/find-active-xrefs-urs-taxids.sql') \
  | find_chunks \
  | splitCsv \
  | combine(data) \
  | export_chunk

  export_chunk.out.count | collect | create_release_note | set { release_note }
  export_chunk.out.xml | collect | set { xml }

  atomic_publish(release_note, xml)
}

workflow {
  search_export()
}
