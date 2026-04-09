#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process create_collection {

  input:
    val(release)

  output:
    val true, emit: collection_ready

  """
  rnac huggingface create-collection $release
  """
}


process create_active_sequences {

  input:
    path(query)

  output:
    path "active_sequences.parquet", emit: parquet

  """
  rnac ftp-export sequences parquet $query active_sequences.parquet
  """
}

process upload_active_sequences {

  input:
    path parquet_file
    val release_num
    val ready_signal

  """
  rnac huggingface create-dataset $release
  rnac huggingface upload-data $release $parquet_file
  """
}

process upload_readme {
  input:
    val release_num
    path update_specs
    val ready_signal

  """
  rnac huggingface create-readme $release $update_specs
  rnac huggingface upload-data $release README.md
  """
}


workflow huggingface {
  take: _flag
  main:
    if (params.export.huggingface.run) {

      update_specs = Channel.fromPath(params.database_update_specs)
      release = Channel.value(params.release)
      active_sql = Channel.fromPath('files/ftp-export/sequences/active.sql')

      create_collection(release)
      collection_ready = create_collection.out.collection_ready

      create_active_sequences(active_sql)
      parquet_ch = create_active_sequences.out.parquet

      upload_active_sequences(parquet_ch, release, collection_ready)
      upload_readme(release, update_specs, collection_ready)

    }
}


workflow {
  huggingface(Channel.of('ready'))
}
