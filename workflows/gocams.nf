#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_urs_taxids {
  input:
  path(query)

  output:
  path('urs_taxids.txt')

  script:
  """
  psql -v ON_ERROR_STOP=1 -At -f "$query" "\$PGDATABASE" > urs_taxids.txt
  """
}

process search_gocams {
  memory params.gocams.memory

  input:
  path(urs_taxids)

  output:
  path("out/go_cam_memberships.${params.writer_format}")

  script:
  """
  rnac gocams search $urs_taxids out
  """
}

process store_gocams {
  input:
  path("go_cam_memberships.${params.writer_format}")
  path(schema)
  path(ctl)
  path(post_load)

  when: params.gocams.load

  script:
  // load-parquet does not create its staging table, so schema.sql always runs.
  // It skips --post-load when nothing loaded, so an empty search never wipes
  // go_cam_memberships; the csv path has no such gate, post-load.sql guards it.
  if (params.writer_format == 'parquet')
    """
    psql -v ON_ERROR_STOP=1 -f $schema "\$PGDATABASE"
    load-parquet load_go_cam_memberships 'go_cam_memberships.parquet' \\
      --post-load $post_load
    """
  else
    """
    psql -v ON_ERROR_STOP=1 -f $schema "\$PGDATABASE"
    cp go_cam_memberships.csv gocams-memberships.csv
    split-and-load $ctl 'gocams-memberships*.csv' ${params.import_data.chunk_size} gocams-memberships
    psql -v ON_ERROR_STOP=1 -f $post_load "\$PGDATABASE"
    """
}

workflow gocams {
  take: _flag
  main:
    if ( !params.gocams.run ) {
      channel.of('gocams skipped') | set { done }
    } else {
      def schema = file('files/gocams/schema.sql')
      def ctl = file('files/gocams/gocams.ctl')
      def post_load = file('files/gocams/post-load.sql')

      channel.fromPath(params.gocams.query) \
      | find_urs_taxids \
      | search_gocams \
      | set { data }

      store_gocams(data, schema, ctl, post_load)
      data | map { _v -> 'gocams done' } | set { done }
    }
  emit: done
}

workflow {
  gocams(channel.of('ready'))
}
