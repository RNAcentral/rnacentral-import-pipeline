process fetch_and_process {
  input:
  path(metadata_query)

  output:
  path('*.{csv,parquet}')

  when: params.databases.ribovision?.run

  script:
  """
  psql -f "$metadata_query" "\$PGDATABASE" > metadata.json
  git clone --depth 1 "$params.databases.ribovision.r2dt_repo" r2dt
  rnac ribovision r2dt-to-fasta r2dt/data sequences.fasta
  rnac ribovision parse metadata.json sequences.fasta
  """
}

workflow ribovision {
  main:
    channel.fromPath('files/import-data/ribovision/metadata.sql') \
    | fetch_and_process \
    | set { data }
  emit: data
}
