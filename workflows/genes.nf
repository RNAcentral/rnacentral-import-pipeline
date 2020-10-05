process get_species {
  input:
  path(query)

  output:
  path('species.csv')

  """
  """
}

process extract_sequences {
  input:
  tuple val(assembly_id), path(query)

  output:
  path('sequences/*.json')

  """
  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id $query $PGDATABASE > raw.json
  rnc genes split raw.json sequences/
  """
}

process build_genes {
  input:
  path(data_file)

  output:
  path('genes.csv')

  """
  rnc genes build $assembly_id $data_file genes.csv
  """
}

process load_data {
  input:
  tuple path('genes*.csv'), path(ctl)

  """
  split-and-load $ctl 'genes*.csv' ${params.import_data.chunk_size} gene-data
  """
}

workflow build_genes {
  get_species.output.set { species }

  species.combine(Channel.fromPath('files/genes/get_species.sql')) |\
  extract_sequences |\
  build_genes |\
  collect() |\
  load_data
}
