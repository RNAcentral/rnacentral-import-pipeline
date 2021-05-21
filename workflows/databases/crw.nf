process fetch_and_process {
  when { params.databases.crw.run }

  input:
  path(metadata_query)

  output:
  path('*.csv')

  """
  psql -f "$metdata_query" "$PGDATABASE" > metadata.json
  git clone "$params.databases.crw.r2dt_repo" r2dt
  rnac crw r2dt-to-fasta r2dt/data/crw-fasta sequences.fasta
  rnac crw parse metadata.json $sequences
  """
}

workflow crw {
  emit: data
  main:
    Channel.fromPath('files/import-data/crw/metadata.sql') \
    | fetch_and_process \
    | set { data }
}
