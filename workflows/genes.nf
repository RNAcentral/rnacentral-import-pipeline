nextflow.enable.dsl=2

process get_species {
  input:
  path(query)

  output:
  path('species.csv')

  """
  psql -v ON_ERROR_STOP=1 $query $PGDATABASE > species.csv
  """
}

process extract_sequences {
  input:
  tuple val(assembly_id), path(query)

  output:
  path('sequences.json')

  """
  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id $query $PGDATABASE > sequences.json
  """
}

process build_genes {
  input:
  path(data_file)

  output:
  path('genes.csv')

  """
  rnc genes build $data_file genes.csv
  """
}

process load_data {
  input:
  tuple path('genes*.csv'), path(ctl), path(post)

  """
  split-and-load $ctl 'genes*.csv' ${params.import_data.chunk_size} gene-data
  psql -v ON_ERROR_STOP=1 -f $post $PGDATABASE
  """
}

workflow build_genes {

  Channel.fromPath('files/genes/species.sql') \
    | get_species \
    | splitCsv \
    | extract_sequences \
    | build_genes \
    | collect \
    | combine(Channel.fromPath('files/genes/load.ctl')) \
    | combine(Channel.fromPath('files/genes/post.ctl')) \
    | load_data
}

workflow {
  build_genes()
}
