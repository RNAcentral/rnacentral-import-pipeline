#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

/* Copy the R2DT models out of the image so r2dt-scan.nf can bind them from the host.
   The image keeps them as data/<source>/cms; r2dt-scan.nf expects cms/<source>. */
process get_r2dt_data {
  container params.r2dt.container
  containerOptions "${params.r2dt_container}"

  input:
    val(data_dir)

  script:
  """
  rm -rf $data_dir/cms
  mkdir -p $data_dir/cms/rfam $data_dir/cms/crw $data_dir/cms/gtrnadb
  cp /rna/r2dt/data/rfam/cms/all.cm $data_dir/cms/rfam/
  cp /rna/r2dt/data/crw/all.cm $data_dir/cms/crw/
  cp /rna/r2dt/data/gtrnadb/cms/*.cm $data_dir/cms/gtrnadb/
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
      def msg = workflow.success
        ? "Environment preparation completed"
        : "Environment preparation failed"
      slack_closure(msg)
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e}"
    }
}
