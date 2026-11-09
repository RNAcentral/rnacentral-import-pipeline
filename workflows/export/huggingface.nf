#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { active_sequences } from './active-sequences'
include { deduplicated_sequences } from './deduplicated'

process create_collection {

  input:
    val(release)

  output:
    val true, emit: collection_ready

  script:
  """
  rnac huggingface create-collection $release
  """
}


process create_parquet {
  tag { name }

  input:
    tuple val(name), path(json)

  output:
    path "${name}.parquet", emit: parquet

  script:
  """
  rnac ftp-export sequences parquet-from-json $json ${name}.parquet
  """
}

process create_dataset {

  input:
    val release_num
    val ready_signal

  output:
    val true, emit: dataset_ready

  script:
  """
  rnac huggingface create-dataset $release_num
  """
}

process upload_active_sequences {

  input:
    path parquet_file
    val release_num
    val ready_signal

  script:
  """
  # container home is read-only; xet writes its cache and logs under HF_HOME
  export HF_HOME=\$PWD/.hf
  rnac huggingface upload-data $release_num $parquet_file
  """
}

process upload_readme {
  input:
    val release_num
    path update_specs
    val ready_signal

  script:
  """
  export HF_HOME=\$PWD/.hf
  rnac huggingface create-readme $release_num $update_specs ${params.clustering.export ? '--with-deduplicated' : ''}
  rnac huggingface upload-data $release_num README.md
  """
}


workflow huggingface {
  take:
    _flag
    active_json
    deduplicated_json
  main:
    if (params.export.huggingface.run) {

      // The README embeds the hand-written per-database version breakdown; without it
      // the dataset card is missing the release context, so fail rather than publish.
      if (!params.export.huggingface.update_specs) {
        error "export.huggingface.update_specs must point to the release notes file"
      }
      update_specs = Channel.fromPath(params.export.huggingface.update_specs, checkIfExists: true)
      release = Channel.value(params.release)

      create_collection(release)
      collection_ready = create_collection.out.collection_ready

      // The deduplicated set is a second split of the same dataset (see create-readme).
      dumps = channel.of('active_sequences').combine(active_json)
        .mix(channel.of('deduplicated').combine(deduplicated_json))
      create_parquet(dumps)
      parquet_ch = create_parquet.out.parquet

      // Dataset repo must exist before either upload; gate both uploads on it so
      // the README upload can't race ahead of dataset creation.
      create_dataset(release, collection_ready)
      dataset_ready = create_dataset.out.dataset_ready.first()

      upload_active_sequences(parquet_ch, release, dataset_ready)
      upload_readme(release, update_specs, dataset_ready)

    }
}


// Escape hatch for re-running HuggingFace on its own. The normal path is export/ftp.nf,
// which shares its dump; running this directly repeats the slow active-sequences query.
workflow {
  active_json = active_sequences()
  deduplicated_sequences(active_json)
  huggingface(Channel.of('ready'), active_json, deduplicated_sequences.out.dedup_json)
}
