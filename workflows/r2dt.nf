nextflow.enable.dsl=2

process import_model_info {
  when { params.r2dt.run }

  output:
  path('info.csv')

  """
  find cms/rfam -type f -name '*.cm' | xargs -I {} cmstat {}  | grep -v ^\\# | awk '{ printf("%s,%d,%d\n", \$3, \$6, \$8); }' info.csv
  """
}

process fetch_model_mapping {
  when { params.r2dt.run }

  input:
  val(_flag)
  path(query)

  output:
  path('mapping.json')

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > mapping.json
  """
}

process extract_sequences {
  when { params.r2dt.run }
  clusterOptions '-sp 100'

  input:
  val(_flag)
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
  rnac r2dt prepare-s3 --allow-missing $mapping $output for-upload file-list
  update-svg.sh file-list $params.r2dt.s3.env
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
  path('data*.csv')
  path(ctl)
  path('attempted*.csv')
  path(attempted_ctl)
  path(urs_sql)
  path(model)
  path(should_show_ctl)

  output:
  val('r2dt done')

  """
  split-and-load $ctl 'data*.csv' ${params.r2dt.data_chunk_size} r2dt-data
  split-and-load $attempted_ctl 'attempted*.csv' ${params.r2dt.data_chunk_size} r2dt-attempted

  psql -f "$urs_sql" "$PGDATABASE" > urs.txt
  rnac r2dt should-show compute $model urs.txt should-show.csv
  split-and-load $should_show_ctl 'should-show*.csv' ${params.r2dt.data_chunk_size} r2dt-should-show
  """
}

workflow common {
  take: ready
  emit: mapping
  main:
    Channel.fromPath('files/r2dt/model_mapping.sql') | set { query }

    fetch_model_mapping(ready, query) | set { mapping }
}

workflow for_database {
  take: sequences
  emit: 
    parsed
    layouts
  main: 
    common | set { model_mapping }

    sequences \
    | flatten \
    | filter_done_md5 \
    | flatten \
    | layout_sequences \
    | combine(model_mapping) \
    | set { layouts }

    layouts | parse_layout | set { parsed }
}

workflow r2dt {
  take: ready
  emit: done
  main:
    Channel.fromPath("files/r2dt/find-sequences.sql") | set { sequences_sql }
    Channel.fromPath('files/r2dt/should-show/model.joblib') | set { ss_model }
    Channel.fromPath('files/r2dt/should-show/query.sql') | set { ss_query }
    Channel.fromPath('files/r2dt/should-show/update.ctl') | set { ss_ctl }
    Channel.fromPath('files/r2dt/load.ctl') | set { load_ctl }
    Channel.fromPath('files/r2dt/attempted.ctl') | set { attempted_ctl }

    ready | common | set { model_mapping }

    extract_sequences(ready, sequences_sql) \
    | flatten \
    | filter { f -> !f.empty() } \
    | split_sequences \
    | flatten \
    | filter { f -> !f.empty() } \
    | layout_sequences \
    | combine(model_mapping) \
    | set { data }

    data | publish_layout
    data | parse_layout 

    parse_layout.out.data | collect | set { data }
    parse_layout.out.attempted | collect | set { attempted }

    store_secondary_structures(data, load_ctl, attempted, attempted_ctl, ss_query, ss_model, ss_ctl) | set { done }
}

workflow {
  r2dt(Channel.from('ready'))
}
