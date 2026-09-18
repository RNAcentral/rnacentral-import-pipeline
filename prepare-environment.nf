#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

/* Copy the R2DT models out of the image so r2dt-scan.nf can bind them from the host */
process get_r2dt_data {
  container params.r2dt.container
  containerOptions "${params.r2dt_container}"

  input:
    val(data_dir)

  script:
  """
  mkdir -p $data_dir
  rm -rf $data_dir/cms
  cp -r /rna/r2dt/data/cms $data_dir/cms
  """
}

workflow prepare_environment {
  main:
    channel.of("Starting environment preparation") | slack_message

    channel.of("$params.r2dt.cms_path/../")| get_r2dt_data
}

workflow {
  main:
    channel.of("Starting...") | slack_message
    prepare_environment()

  onComplete:
    try {
      slack_closure("Environment preparation completed")
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e}"
    }
}
