#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process parse_clustering {
  input:
    path(parquet_dir)

  output:
    path('parsed')

  script:
  """
  rnac clustering parse $parquet_dir parsed
  """
}

process store_clustering {
  memory 16.GB

  input:
    path(parsed)
    path(schema)
    path(post_load)

  script:
  """
  psql -v ON_ERROR_STOP=1 -f $schema "\$PGDATABASE"
  load-parquet load_sequence_cluster '$parsed/clusters.parquet' --truncate
  load-parquet load_sequence_absorption '$parsed/absorption.parquet' --truncate --post-load $post_load
  """
}

workflow clustering {
  take: _flag
  main:
    if (params.clustering.load) {
      if (!params.clustering.dir) {
        error "clustering.load is set but clustering.dir is not"
      }
      parse_clustering(file(params.clustering.dir, checkIfExists: true))
      store_clustering(
        parse_clustering.out,
        file('files/clustering/schema.sql'),
        file('files/clustering/post-load.sql'),
      )
    }
}

workflow {
  clustering(channel.of('ready'))
}
