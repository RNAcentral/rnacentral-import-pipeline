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

process extract_data {
  tag { "$assembly_id" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  maxForks params.genes.extract_sequences.maxForks

  input:
  tuple val(assembly_id), val(taxid), path(query), path(counts_query), path(genes_query)

  output:
  tuple val(assembly_id), path('sequences.json'), path('counts.json'), path("pseudo.json"), path('repetitive.bed')

  """
  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id -v taxid=$taxid -f $query $PGDATABASE > sequences.json
  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id -v taxid=$taxid -f $counts_query $PGDATABASE > counts.json
  psql -v ON_ERROR_STOP=1 -v assembly_id=$assembly_id -f $genes_query "$PGDATABASE" > pseudo.json
  touch repetitive.bed
  """
}

process build {
  tag { "$assembly_id" }
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"
  errorStrategy 'ignore'

  input:
  tuple val(assembly_id), path(data_file), path(counts), path(genes), path(repetative)

  output:
  path 'locus.csv', emit: locus
  path 'rejected.csv', emit: rejected

  """
  rnac genes build $data_file $counts $genes $repetative .
  """
}

process load_data {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path('locus*.csv')
  path('rejected*.csv')
  path(load_locus)
  path(load_status)
  path(post_load)

  """
  split-and-load $load_locus 'locus*.csv' ${params.import_data.chunk_size} locus-data
  STATUS='rejected' split-and-load $load_status 'rejected*.csv' ${params.import_data.chunk_size} rejected-data
  psql -v ON_ERROR_STOP=1 -f $post_load $PGDATABASE
  """
}

workflow genes {
  Channel.fromPath('files/genes/load-status.ctl').set { load_status }
  Channel.fromPath('files/genes/load.ctl').set { load }
  Channel.fromPath('files/genes/post-load.sql').set { post_load }

  Channel.fromPath('files/genes/species.sql') \
    | combine(Channel.fromPath('files/genes/schema.sql')) \
    | get_species \
    | splitCsv \
    | combine(Channel.fromPath('files/genes/data.sql')) \
    | combine(Channel.fromPath('files/genes/counts.sql')) \
    | combine(Channel.fromPath('files/genes/pseudogenes.sql')) \
    | extract_data \
    | build

  build.out.locus | collect | set { locus }
  build.out.rejected | collect | set { rejected }

  load_data(locus, rejected, load, load_status, post_load)
}

workflow {
  genes()
}
