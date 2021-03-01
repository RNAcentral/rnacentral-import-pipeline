process find_chunks {
  executor 'local'
  when: params.ftp_export.ensembl.run

  output:
  path('chunks')

  """
  rnac upi-ranges ${params.ftp_export.ensembl.chunk_size} > chunks
  """
}

process query_chunk {
  maxForks params.ftp_export.ensembl.maxForks

  input:
  tuple val(table), val(min), val(max), path(query)

  output:
  tuple val(min), val(max), path('raw_xrefs.json')

  """
  psql -f $query --variable min=$min --variable max=$max "$PGDATABASE" > raw_xrefs.json
  """
}

process process_chunk {
  publishDir "${params.ftp_export.publish}/json/", mode: 'copy'

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
  | combine(schema) \
  | process_chunk
}
