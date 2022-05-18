#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

/* Get some data downloaded and in the right place */

process wget_r2dt_data{
  errorStrategy { slack_closure("""Wget is having trouble downloading r2dt data.
  This process will retry a maximum of 5 times;
  if the problem persists, try rerunning with the
  --use_datamover option""".stripIndent());
  return 'retry' }
  maxRetries 5

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
  wget --read-timeout=30 -t 1 -O cms.tar.gz http://ftp.ebi.ac.uk/pub/databases/RNAcentral/r2dt/1.2/cms.tar.gz

  tar -xf cms.tar.gz

  """
}


/* On the cluster this is much much faster than wget */
process dmget_r2dt_data {
  queue 'datamover'
  executor 'lsf'
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

  cp /nfs/ftp/public/databases/RNAcentral/r2dt/1.2/cms.tar.gz .

  tar -xf cms.tar.gz
  """
}

/*
This is a bit opaque. It creates two channels, one of which is empty depending
on the value of params.use_datamover. Processes consuming empty channels are not
executed.

I *think* this is preferable to using a when directive in the process, but not
100% sure
*/

workflow prepare_environment {
  main:
    Channel.of("Starting environment preparation") | slack_message
    (dm_inch, wg_inch) = ( params.use_datamover ?
            [Channel.of("$workflow.launchDir/singularity/bind/r2dt/data"), Channel.empty()]
          : [Channel.empty(), Channel.of("$workflow.launchDir/singularity/bind/r2dt/data")] )

    dm_inch | dmget_r2dt_data
    wg_inch | wget_r2dt_data

}

workflow {
  prepare_environment()
}

workflow.onComplete {
  slack_closure("Environment preparation completed")
}
