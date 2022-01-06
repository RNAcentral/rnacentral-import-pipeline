process find_chunks {
  executor 'local'

  output:
  path('chunks')

  """
  rnac upi-ranges ${params.export.ftp.ensembl.chunk_size} > chunks
  """
}

process query_chunk {
  maxForks params.export.ftp.ensembl.maxForks

  input:
  tuple val(table), val(min), val(max), path(query)

  output:
  tuple val(min), val(max), path('raw_xrefs.json')

  """
  psql -f $query --variable min=$min --variable max=$max "$PGDATABASE" > raw_xrefs.json
  """
}

process process_chunk {
  publishDir "${params.export.ftp.publish}/json/", mode: 'copy'

  input:
  tuple val(min), val(max), file(raw), path(schema)

  output:
  path("ensembl-xref-$min-${max}.json")

  """
  rnac ftp-export ensembl --schema=$schema $raw ensembl-xref-$min-${max}.json
  """
}

workflow ensembl_export {
  Channel.fromPath('files/ftp-export/ensembl/schema.json') | set { schema }
  Channel.fromPath('files/ftp-export/ensembl/ensembl-xrefs.sql') | set { query }

  find_chunks \
  | splitCsv \
  | combine(query) \
  | query_chunk \
  | filter { _min, _max, xrefs -> !xrefs.isEmpty() } \
  | combine(schema) \
  | process_chunk
}
