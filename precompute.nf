#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_genomes_with_repeats {
  input:
  path(query)

  output:
  path("assemblies.csv")

  """
  psql -f $query "$PGDATABASE" > assemblies.csv
  """
}

process fetch_ensembl_data {
  input:
  tuple val(species), val(assembly), val(host)

  output:
  path("repeat-${assembly}.pickle"), emit: repeats

  """
  rnac repeats url-for $species $assembly $host - | xargs -I {} fetch generic '{}' ${species}.fasta.gz
  gzip -d ${species}.fasta.gz
  rnac repeats compute-ranges $assembly ${species}.fasta repeat-${assembly}.pickle
  """
}

process build_precompute_context {
  input:
  path('repeats*.pickle')

  output:
  path('context.pickle')

  """
  rnac repeats build-tree repeats*.pickle repeat-tree.pickle
  rnac precompute build-context repeat-tree.pickle pseudo-tree.pickle context.pickle
  """
}

process find_precompute_upis {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  path(sql)

  output:
  path('ranges.txt')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=$params.precompute.tablename \
    -f "$sql" "$PGDATABASE"
  rnac upi-ranges --table-name $params.precompute.tablename ${params.precompute.max_entries} ranges.txt
  """
}

process precompute_range_query {
  tag { "$min-$max" }
  beforeScript 'slack db-work precompute-range || true'
  afterScript 'slack db-done precompute-range || true'
  maxForks params.precompute.maxForks
  container ''

  input:
  tuple val(tablename), val(min), val(max), path(query)

  output:
  tuple val(min), val(max), path('raw-precompute.json')

  """
  psql \
    --variable ON_ERROR_STOP=1 \
    --variable tablename=${params.precompute.tablename} \
    --variable min=$min \
    --variable max=$max \
    -f "$query" \
    '$PGDATABASE' > raw-precompute.json
  """
}

process process_range {
  tag { "$min-$max" }
  memory params.precompute.range.memory
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  input:
  tuple val(min), val(max), path(raw), path(context)

  output:
  path('precompute.csv'), emit: data
  path('qa.csv'), emit: qa

  """
  rnac precompute from-file $context $raw
  """
}

process load_precomputed_data {
  beforeScript 'slack db-work loading-precompute || true'
  afterScript 'slack db-done loading-precompute || true'
  container ''

  input:
  tuple path('precompute*.csv'), path('qa*.csv'), path(pre_ctl), path(qa_ctl), path(post)

  script:
  def tablename = params.precompute.tablename
  """
  split-and-load $pre_ctl 'precompute*.csv' ${params.import_data.chunk_size} precompute
  split-and-load $qa_ctl 'qa*.csv' ${params.import_data.chunk_size} qa
  psql -v ON_ERROR_STOP=1 -f $post "$PGDATABASE"
  psql -v ON_ERROR_STOP=1 -c 'DROP TABLE IF EXISTS $tablename' "$PGDATABASE"
  """
}

workflow precompute {
  take:
    path(method)

  main:
    find_genomes_with_repeats(Channel.fromPath('files/repeat/find-assembiles.sql')) \
    | splitCsv() \
    | fetch_repeats

    fetch_ensembl_data.repeats | collect | set { repeats }

    build_precompute_context(repeats) | set { context }

    find_ranges(method) \
    | splitCsv() \
    | combine(Channel.fromPath('files/precompute/query.sql')) \
    | query_range \
    | combine(context) \
    | process_range \

    process_range.out.data | collect | set { data }
    process_range.out.qa | collect | set { qa }

    data \
    | combine(qa) \
    | combine(Channel.fromPath('files/precompute/load.ctl')) \
    | combine(Channel.fromPath('files/precompute/qa.ctl')) \
    | combine(Channel.fromPath('files/precompute/post-load.sql')) \
    | load_data
}

workflow {
  precompute(Channel.fromPath("files/precompute/methods/${params.precompute.method.replace('_', '-')}.sql"))
}
