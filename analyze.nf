#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { cpat } from './workflows/cpat'
include { genome_mapping } from './workflows/genome-mapping'
include { r2dt } from './workflows/r2dt'
include { rfam_scan } from './workflows/rfam-scan'
include { rediportal } from './workflows/databases/rediportal'

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

workflow analyze {
  take: ready
  emit: done
  main:
    Channel.of("Starting analyze pipeline") | slack_message
    ready | (genome_mapping & rfam_scan & r2dt & cpat ) | mix | collect | rediportal | set { done }
}

workflow {
  analyze(Channel.of('ready'))
}


workflow.onComplete {
  slack_closure("Analyze workflow completed")

}

workflow.onError {

  slack_closure("Analyze workflow hit an error and crashed")

}
