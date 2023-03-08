#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

/* Get some data downloaded and in the right place */

/* On the cluster this is much much faster than wget */
process get_r2dt_data {
  queue 'datamover'
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

  cp /nfs/ftp/public/databases/RNAcentral/r2dt/1.3/cms.tar.gz .

  tar -xf cms.tar.gz --strip-components=1 -C ./cms
  """
}

workflow prepare_environment {
  main:
    Channel.of("Starting environment preparation") | slack_message

    Channel.of("$params.r2dt.cms_path/../")| get_r2dt_data
}

workflow {
  Channel.of("Starting...") | slack_message
  prepare_environment()
}

workflow.onComplete {
  slack_closure("Environment preparation completed")
}
