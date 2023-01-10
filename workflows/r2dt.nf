nextflow.enable.dsl=2

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
    -v 'sequence_count=${params.r2dt.sequence_count}' \
    -f "$query" "$PGDATABASE" > raw.json
  mkdir parts/
  split --number=l/4000 --additional-suffix='.json' raw.json parts/
  """
}

process split_sequences {

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
  errorStrategy { task.exitStatus = 130 ? 'ignore' : 'finish' }

  input:
  path(sequences)

  output:
  tuple path("$sequences"), path('output'), path('version')

  """
  esl-sfetch --index $sequences
  r2dt.py draw $sequences output/
  r2dt.py version | perl -ne 'm/(\\d\\.\\d)/ && print "\$1\\n"' > version
  """
}

process publish_layout {
  maxForks 50
  errorStrategy { task.attempt < 5 ? "retry" : "ignore" }
  maxRetries 5
  queue 'datamover'

  input:
  tuple path(sequences), path(output), path(_version), path(mapping)

  output:
  val 'done', emit: flag

  """
  rnac r2dt publish --allow-missing $mapping $output $params.r2dt.publish
  rnac r2dt prepare-s3 --allow-missing $mapping $output for-upload file-list
  update-svg.sh file-list $params.r2dt.s3.env
  """
}

process parse_layout {
  input:
  tuple path(sequences), path(to_parse), path(version), path(mapping)
  errorStrategy "ignore"

  output:
  path "data.csv", emit: data
  path 'attempted.csv', emit: attempted

  """
  rnac r2dt process-svgs --allow-missing $mapping $to_parse data.csv
  rnac r2dt create-attempted $sequences $version attempted.csv
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
  val(_flag)

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

    publish_layout.out.flag | collect | map { _ -> 'ready' } | set { uploaded }

    store_secondary_structures(data, load_ctl, attempted, attempted_ctl, ss_query, ss_model, ss_ctl, uploaded) | set { done }
}

workflow {
  r2dt(Channel.from('ready'))
}
