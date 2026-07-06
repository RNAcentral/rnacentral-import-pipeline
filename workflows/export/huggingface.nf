#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { active_sequences } from './active-sequences'

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
    path(active_json)

  output:
    path "active_sequences.parquet", emit: parquet

  """
  rnac ftp-export sequences parquet-from-json $active_json active_sequences.parquet
  """
}

process create_dataset {

  input:
    val release_num
    val ready_signal

  output:
    val true, emit: dataset_ready

  """
  rnac huggingface create-dataset $release_num
  """
}

process upload_active_sequences {

  input:
    path parquet_file
    val release_num
    val ready_signal

  """
  rnac huggingface upload-data $release_num $parquet_file
  """
}

process upload_readme {
  input:
    val release_num
    path update_specs
    val ready_signal

  """
  rnac huggingface create-readme $release_num $update_specs
  rnac huggingface upload-data $release_num README.md
  """
}


workflow huggingface {
  take:
    _flag
    active_json
  main:
    if (params.export.huggingface.run) {

      update_specs = Channel.fromPath(params.export.huggingface.update_specs)
      release = Channel.value(params.release)

      create_collection(release)
      collection_ready = create_collection.out.collection_ready

      create_active_sequences(active_json)
      parquet_ch = create_active_sequences.out.parquet

      // Dataset repo must exist before either upload; gate both uploads on it so
      // the README upload can't race ahead of dataset creation.
      create_dataset(release, collection_ready)
      dataset_ready = create_dataset.out.dataset_ready.first()

      upload_active_sequences(parquet_ch, release, dataset_ready)
      upload_readme(release, update_specs, dataset_ready)

    }
}


workflow {
  huggingface(Channel.of('ready'), active_sequences())
}
