nextflow.enable.dsl=2

include { model_info } from './r2dt/model-info.nf'

process fetch_model_mapping {
  memory '16 MB'

  input:
  val(_flag)
  path(query)

  output:
  path('mapping.json')

  when: params.r2dt?.run

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $query "\$PGDATABASE" > mapping.json
  """
}

process get_partitions {
  input:
  val(_flag)
  path(query)

  output:
  path('databases.csv')

  when: params.r2dt?.run

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -f $query "\$PGDATABASE" \
    > databases.csv
  """
}


process fetch_xrefs {
  input:
  tuple val(partition), path(query)

  output:
  path('urs_xref.csv')

  when: params.r2dt?.run

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -v tablename=xref_p${partition}_not_deleted \
    -f $query "\$PGDATABASE" \
    > urs_xref.csv
  """
}

process fetch_tracked {
  input:
  val(_flag)
  path(query)

  output:
  path('urs_tracked.csv')

  when: params.r2dt?.run

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -f $query "\$PGDATABASE" \
    > urs_tracked.csv
  """
}


process extract_sequences {
  memory '12GB'

  input:
  path(xrefs)
  path(tracked)
  path(query)

  output:
  path('raw.json')

  when: params.r2dt?.run

  script:
  """
    rnac r2dt prepare-sequences $xrefs $tracked urs_to_fetch.csv

    psql \
      -v ON_ERROR_STOP=1 \
      -v max_len=10000 \
      -v 'sequence_count=${params.r2dt.sequence_count}' \
      -q \
      -f $query \$PGDATABASE > raw.json
  """
}

process split_sequences {

  memory { 8.GB * task.attempt }
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
  maxForks params.r2dt.layout_max_forks
  memory params.r2dt.layout.memory
  container params.r2dt.container
  containerOptions "${params.r2dt_container}"
  errorStrategy 'ignore'

  input:
  path(sequences)

  output:
  tuple path("$sequences"), path('output'), path('version')

  script:
  """
  esl-sfetch --index $sequences
  r2dt.py draw $sequences output/ || true
  r2dt.py version | perl -ne 'm/(\\d\\.\\d)/ && print "\$1\\n"' > version
  """
}

process publish_layout {
  maxForks params.r2dt.publish_max_forks
  errorStrategy { task.attempt < 5 ? "retry" : "finish" }
  maxRetries 5
  queue 'datamover'
  memory { 1.GB * task.attempt }

  input:
  tuple path(_sequences, stageAs: 'seq_*.fasta'), path(outputs, stageAs: 'output_*'), path(_versions, stageAs: 'version_*'), path(mapping)

  output:
  val 'done', emit: flag
  // Only written when an upload failed and was tolerated, so the channel stays
  // empty on a clean run.
  path 'upload-failures.txt', optional: true, emit: failures

  script:
  // prepare-s3 opens its file list with 'w', so each chunk needs its own to concatenate.
  """
  outs=( $outputs )

  for i in "\${!outs[@]}"; do
    rnac r2dt publish --allow-missing $mapping "\${outs[\$i]}" $params.r2dt.publish
    rnac r2dt prepare-s3 --allow-missing $mapping "\${outs[\$i]}" for-upload "file-list.\$i"
  done

  cat file-list.* > file-list
  rnac r2dt upload-s3 --env $params.r2dt.s3.env \\
    --allow-failures ${params.r2dt.allow_upload_failures} \\
    --failure-list upload-failures.txt \\
    file-list
  """
}

process parse_layout {
  maxForks params.r2dt.parse_max_forks
  memory '2 GB'
  errorStrategy "ignore"

  input:
  tuple path(sequences, stageAs: 'seq_*.fasta'), path(to_parse, stageAs: 'output_*'), path(versions, stageAs: 'version_*'), path(mapping)

  output:
  path "data_*.csv", emit: data, optional: true
  path "attempted_*.${params.writer_format}", emit: attempted, optional: true

  script:
  // Skip a failed chunk rather than the batch, so errorStrategy 'ignore' still costs one chunk.
  """
  seqs=( $sequences )
  outs=( $to_parse )
  vers=( $versions )

  for i in "\${!outs[@]}"; do
    if ! rnac r2dt process-svgs --allow-missing $mapping "\${outs[\$i]}" "data_\$i.csv"; then
      echo "process-svgs failed for \${outs[\$i]}; skipping chunk" >&2
      rm -f "data_\$i.csv"
      continue
    fi
    rnac r2dt create-attempted "\${seqs[\$i]}" "\${vers[\$i]}" "attempted_\$i.${params.writer_format}"
  done
  """
}

process store_secondary_structures {
  memory 9.GB

  input:
  path('data*.csv')
  path(ctl)
  path("attempted*.${params.writer_format}")
  path(attempted_ctl)
  path(attempted_post_load)
  path(urs_sql)
  path(model)
  path(should_show_ctl)
  path(upload_failures)
  val(_flag)

  output:
  val('r2dt done')

  script:
  // Hits side (data*.csv) and should-show stay on the legacy pgloader path
  // for now; only the attempted side has been migrated to parquet.
  def attempted_cmd = (params.writer_format == 'parquet') ?
    """
    load-parquet load_traveler_attempted 'attempted*.parquet' \\
      --truncate \\
      --post-load $attempted_post_load
    """ :
    "split-and-load $attempted_ctl 'attempted*.csv' ${params.r2dt.data_chunk_size} r2dt-attempted"
  """
  rnac r2dt drop-failed-uploads $upload_failures \\
    data*.csv attempted*.${params.writer_format}

  split-and-load $ctl 'data*.csv' ${params.r2dt.data_chunk_size} r2dt-data
  $attempted_cmd

  psql -f "$urs_sql" "\$PGDATABASE" > urs.txt
  rnac r2dt should-show compute $model urs.txt should-show.csv
  split-and-load $should_show_ctl 'should-show*.csv' ${params.r2dt.data_chunk_size} r2dt-should-show
  """
}

workflow common {
  take: ready
  main:
    channel.fromPath('files/r2dt/model_mapping.sql') | set { query }

    fetch_model_mapping(ready, query) | set { mapping }
  emit: mapping
}

workflow r2dt {
  take: ready
  main:
    if (params.r2dt.run) {
      channel.fromPath("files/r2dt/find-sequences.sql") | set { sequences_sql }
      channel.fromPath("files/r2dt/fetch-partition-xref.sql") | set { xref_sql }
      channel.fromPath("files/r2dt/fetch-tracked.sql") | set { tracked_sql }
      channel.fromPath("files/r2dt/fetch-partitions.sql") | set { partitions_sql }
      channel.fromPath('files/r2dt/should-show/model.joblib') | set { ss_model }
      channel.fromPath('files/r2dt/should-show/query.sql') | set { ss_query }
      channel.fromPath('files/r2dt/should-show/update.ctl') | set { ss_ctl }
      channel.fromPath('files/r2dt/load.ctl') | set { load_ctl }
      channel.fromPath('files/r2dt/attempted.ctl') | set { attempted_ctl }
      channel.fromPath('files/r2dt/attempted-post-load.sql') | set { attempted_post_load }

      model_info(ready) | set { models_ready }

      fetch_tracked(models_ready, tracked_sql) | set { tracked }

      get_partitions(models_ready, partitions_sql) \
      | splitCsv \
      | combine(xref_sql) \
      | fetch_xrefs \
      | collectFile \
      | set { partitions }

      models_ready | common | set { model_mapping }

      extract_sequences(partitions, tracked, sequences_sql) \
      | flatten \
      | filter { f -> !f.empty() } \
      | split_sequences \
      | flatten \
      | filter { f -> !f.empty() } \
      | layout_sequences \
      | collate(params.r2dt.batch_size) \
      | map { batch -> tuple(
          batch.collect { row -> row[0] },
          batch.collect { row -> row[1] },
          batch.collect { row -> row[2] },
        ) } \
      | combine(model_mapping) \
      | set { data }

      data | publish_layout
      data | parse_layout

      parse_layout.out.data | collect | set { data }
      parse_layout.out.attempted | collect | set { attempted }

      publish_layout.out.flag | collect | map { _flags -> 'ready' } | set { uploaded }

      // Every URS that could not be uploaded, in one file next to the SVGs.
      // store_secondary_structures drops those rows before loading, so a
      // tolerated upload failure leaves the sequence undone rather than
      // recorded with an SVG that is not in S3. The placeholder keeps the
      // input satisfied on a clean run, where the channel is empty.
      publish_layout.out.failures \
      | collectFile(name: 'upload-failures.txt', storeDir: params.r2dt.publish) \
      | ifEmpty(file('files/r2dt/no-upload-failures.txt')) \
      | set { upload_failures }

      store_secondary_structures(data, load_ctl, attempted, attempted_ctl, attempted_post_load, ss_query, ss_model, ss_ctl, upload_failures, uploaded) | set { done }
    } else {
      channel.of('r2dt skipped') | set { done }
    }
  emit: done
}

workflow {
  r2dt(channel.from('ready'))
}
