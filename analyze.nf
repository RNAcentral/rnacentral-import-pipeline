#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { cpat } from './workflows/cpat'
include { genome_mapping } from './workflows/genome-mapping'
include { r2dt } from './workflows/r2dt.nf'
include { rfam_scan } from './workflows/rfam-scan'
include { rediportal } from './workflows/databases/rediportal'
include { stopfree } from './workflows/stopfree'
include { tcode } from './workflows/tcode'

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

workflow analyze {
  take: ready
  main:
    Channel.of("Starting analyze pipeline") | slack_message
    ready | (genome_mapping & rfam_scan & r2dt & cpat & tcode & stopfree ) | mix | collect | rediportal | set { done }
  emit: done
}

workflow {
  analyze(Channel.of('ready'))

  workflow.onComplete {
    try {
      def msg = workflow?.success ? "Analyze workflow completed" : "Analyze workflow completed with errors"
      slack_closure(msg)
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e.message}"
    }
  }

  workflow.onError {
    try {
      slack_closure("Analyze workflow hit an error and crashed")
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e.message}"
    }
  }
}
