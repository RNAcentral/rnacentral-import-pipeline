nextflow.enable.dsl=2

process fetch_model_mapping {
  when { params.r2dt.run }

  memory '16 MB'

  input:
  val(_flag)
  path(query)

  output:
  path('mapping.json')

  """
  psql -v ON_ERROR_STOP=1 -f $query "$PGDATABASE" > mapping.json
  """
}

process get_partitions {

  input:
  val(_flag)
  path(query)

  output:
  path('databases.csv')

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -f $query "$PGDATABASE" \
    > databases.csv
  """
}


process fetch_xrefs {
  // when { params.r2dt.run }

  input:
  tuple val(partition), path(query)

  output:
  path('urs_xref.csv')

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=xref_p${partition}_not_deleted \
    -f $query "$PGDATABASE" \
    > urs_xref.csv
  """
}

process fetch_tracked {
  // when { params.r2dt.run }

  input:
  tuple val(_flag)
  path(query)

  output:
  path('urs_tracked.csv')

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -f $query "$PGDATABASE" \
    > urs_tracked.csv
  """
}


process extract_sequences {
  // when { params.r2dt.run }

  memory '12GB'

  input:
  path(xrefs)
  path(tracked)
  path(query)

  output:
  path('raw.json')

  script:
  """
    rnac r2dt prepare-sequences $xrefs $tracked urs_to_fetch.csv

    psql \
      -v ON_ERROR_STOP=1 \
      -v max_len=10000 \
      -v 'sequence_count=${params.r2dt.sequence_count}' \
      -q \
      -f $query $PGDATABASE > raw.json
  """
}

process split_sequences {

  memory { 8.MB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

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
  memory '256 MB'

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
    memory '256 MB'
    errorStrategy "ignore"

  input:
  tuple path(sequences), path(to_parse), path(version), path(mapping)

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
    Channel.fromPath("files/r2dt/fetch-partition-xref.sql") | set { xref_sql }
    Channel.fromPath("files/r2dt/fetch-tracked.sql") | set { tracked_sql }
    Channel.fromPath("files/r2dt/fetch-partitions.sql") | set { partitions_sql }
    Channel.fromPath('files/r2dt/should-show/model.joblib') | set { ss_model }
    Channel.fromPath('files/r2dt/should-show/query.sql') | set { ss_query }
    Channel.fromPath('files/r2dt/should-show/update.ctl') | set { ss_ctl }
    Channel.fromPath('files/r2dt/load.ctl') | set { load_ctl }
    Channel.fromPath('files/r2dt/attempted.ctl') | set { attempted_ctl }

    fetch_tracked(ready, tracked_sql) | set { tracked }


    get_partitions(ready, partitions_sql) \
    | splitCsv \
    | combine(xref_sql) \
    | fetch_xrefs \
    | collectFile \
    | set { partitions }

    ready | common | set { model_mapping }

  extract_sequences(partitions, tracked, sequences_sql) \
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
