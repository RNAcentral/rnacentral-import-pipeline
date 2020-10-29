nextflow.enable.dsl=2

process get_species {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple path(query), path(setup)

  output:
  path('species.csv')

  """
  psql -v ON_ERROR_STOP=1 -f $setup $PGDATABASE
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
  path('locus.csv'), emit: locus
  path('rejected.csv'), emit: rejected
  path('ignored.csv'), emit: ignored

  """
  rnac genes build $data_file .
  """
}

process load_data {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('locus*.csv')
  path('rejected*.csv')
  path('ignored*.csv')
  path(load_locus)
  path(load_status)
  path(post_load)

  """
  split-and-load $load_locus 'locus*.csv' ${params.import_data.chunk_size} locus-data
  STATUS='rejected' split-and-load $load_status 'rejected*.csv' ${params.import_data.chunk_size} rejected-data
  STATUS='ignored' split-and-load $load_status 'ignored*.csv' ${params.import_data.chunk_size} ignored-data
  psql -v ON_ERROR_STOP=1 -f $post_load $PGDATABASE
  """
}

workflow build_genes {
  Channel.fromPath('files/genes/load-status.ctl').set { load_status }
  Channel.fromPath('files/genes/load.ctl').set { load }
  Channel.fromPath('files/genes/post-load.sql').set { post_load }

  Channel.fromPath('files/genes/species.sql') \
    | combine(Channel.fromPath('files/genes/schema.sql')) \
    | get_species \
    | splitCsv \
    | map { row -> row[0] } \
    | combine(Channel.fromPath('files/genes/data.sql')) \
    | extract_sequences \
    | build

    build.out.locus | collect | set { locus }
    build.out.rejected | collect | set { rejected }
    build.out.ignored | collect | set { ignored }
    
    load_data(locus, rejected, ignored, load, load_status, post_load)
}

workflow {
  build_genes()
}
