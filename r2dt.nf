nextflow.enable.dsl=2

process fetch_model_mapping {
  input:
  path(query)

  output:
  path('mapping.json')

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > mapping.json
  """
}

process extract_sequences {
  clusterOptions '-sp 100'

  input:
  path(query)

  output:
  path('parts/*.json')

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v 'tablename=${params.r2dt.tablename}' \
    -v max_len=10000 \
    -f "$query" "$PGDATABASE" > raw.json
  mkdir parts/
  split --number=l/4000 --additional-suffix='.json' raw.json parts/
  """
}

process split_sequences {
  clusterOptions '-sp 100'

  input:
  path("raw.json")

  output:
  path('parts/*.fasta')

  script:
  def chunk_size = params.r2dt.sequence_chunk_size
  """
  mkdir parts/
  split --lines=${chunk_size} --additional-suffix='.fasta' --filter 'json2fasta - - >> \$FILE' raw.json parts/
  """
}

process layout_sequences {
  tag { "${sequences}" }
  memory params.r2dt.layout.memory
  container params.r2dt.container
  containerOptions "--bind ${params.r2dt.cms_path}:/rna/r2dt/data/cms" 
  errorStrategy { task.exitStatus = 130 ? 'ignore' : 'terminate' }

  input:
  path(sequences)

  output:
  tuple path("$sequences"), path('output')

  """
  esl-sfetch --index $sequences
  r2dt.py draw $sequences output/
  """
}

process publish_layout {
  input:
  tuple path(sequences), path(output), path(mapping)

  """
  rnac r2dt publish --allow-missing $mapping $output $params.r2dt.publish
  """
}

process parse_layout {
  input:
  tuple path(sequences), path(to_parse), path(mapping)

  output:
  path "data.csv", emit: data
  path 'attempted.csv', emit: attempted

  """
  rnac r2dt process-svgs --allow-missing $mapping $to_parse data.csv
  rnac r2dt create-attempted $sequences attempted.csv
  """
}

process store_secondary_structures {
  memory params.r2dt.store.memory

  input:
  tuple path('data*.csv'), path(ctl)
  tuple path('attempted*.csv'), path(attempted_ctl)

  """
  split-and-load $ctl 'data*.csv' ${params.r2dt.data_chunk_size} traveler-data
  split-and-load $attempted_ctl 'attemped*.csv' ${params.r2dt.data_chunk_size} traveler-attempted
  """
}

workflow common {
  emit: fasta
  main:
    Channel.fromPath('files/r2dt/model_mapping.sql') \
    | fetch_model_mapping \
    | set { mapping }
}

workflow for_database {
  take: sequences
  emit: data
  main: 
    common | set { model_mapping }

    sequences \
    | flatten \
    | filter_done_md5 \
    | flatten \
    | layout_sequences \
    | combine(model_mapping) \
    | set { data }
}

workflow r2dt {
  common | set { model_mapping }

  Channel.fromPath("files/r2dt/find-sequences.sql") \
  | extract_sequences \
  | flatten \
  | split_sequences \
  | flatten \
  | layout_sequences \
  | combine(model_mapping) \
  | set { data }

  data | publish_layout
  data | parse_layout 
  parse_layout.out.data | collect | combine(Channel.fromPath('files/r2dt/load.ctl')) | set { to_load }
  parse_layout.out.attempted | collect | combine(Channel.fromPath('files/r2dt/attempted.ctl')) | set { attempted }

  store_secondary_structures(to_load, attempted)
}

workflow {
  r2dt()
}
