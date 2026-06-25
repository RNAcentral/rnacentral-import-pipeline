#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

/* Get some data downloaded and in the right place */

/* On the cluster this is much much faster than wget */
process get_r2dt_data {
  container ''

  input:
    val(data_dir)

  script:
  """
  echo "$data_dir"
  if [ ! -d $data_dir ]
  then
    mkdir -p $data_dir
  fi

  cd $data_dir

  wget https://github.com/r2dt-bio/R2DT/releases/download/v2.0/cms.tar.gz

  tar -xf cms.tar.gz --strip-components=1 -C ./cms
  """
}

process get_dfam_data {
  container ''

  input:
    tuple val(data_dir), val(url)

  script:
  """
  mkdir -p $data_dir
  cd $data_dir
  fname=\$(basename "$url")
  wget -O "\$fname" "$url"
  gunzip -f "\$fname"
  """
}

workflow prepare_environment {
  main:
    Channel.of("Starting environment preparation") | slack_message

    Channel.of("$params.r2dt.cms_path/../")| get_r2dt_data
    Channel.of([params.repeatmasker.dfam_path, params.repeatmasker.dfam_url]) | get_dfam_data
}

workflow {
  Channel.of("Starting...") | slack_message
  prepare_environment()

  workflow.onComplete {
    try {
      slack_closure("Environment preparation completed")
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e.message}"
    }
  }
}
