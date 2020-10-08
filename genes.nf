nextflow.enable.dsl=2

process get_species {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(query)

  output:
  path('species.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $query $PGDATABASE > species.csv
  """
}

process extract_sequences {
  tag { "$assembly_id" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  maxForks params.genes.extract_sequences.maxForks

  input:
  tuple val(assembly_id), path(query)

  output:
  tuple val(assembly_id), path('sequences.json')

  """
  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id -f $query $PGDATABASE > sequences.json
  """
}

process build {
  tag { "$assembly_id" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(assembly_id), path(data_file)

  output:
  path('genes.csv')

  """
  rnac genes build $data_file genes.csv
  """
}

process load_data {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

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
    | map { row -> row[0] } \
    | combine(Channel.fromPath('files/genes/data.sql')) \
    | extract_sequences \
    | build \
    | collect \
    | map { files -> [files, "$baseDir/files/genes/load.ctl", "$baseDir/files/genes/post-load.sql"] } \
    | load_data
}

workflow {
  build_genes()
}
