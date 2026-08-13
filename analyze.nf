#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { cpat } from './workflows/cpat'
include { genome_mapping } from './workflows/genome-mapping'
include { r2dt } from './workflows/r2dt.nf'
include { rfam_scan } from './workflows/rfam-scan'
include { stopfree } from './workflows/stopfree'
include { tcode } from './workflows/tcode'

include { qc_analyze } from './workflows/utils/qc'

include { slack_closure } from './workflows/utils/slack'
include { slack_message } from './workflows/utils/slack'

workflow analyze {
  take: ready
  main:
    channel.of("Starting analyze pipeline") | slack_message

    // Value channel so start triggers both the analyses and the QC before-snapshot.
    ready | first | set { start }

    start | (genome_mapping & rfam_scan & r2dt & cpat & tcode & stopfree )

    // QC: snapshot counts at start, report per-analysis once all finish.
    channel.empty() \
    | mix(
      genome_mapping.out,
      rfam_scan.out,
      r2dt.out,
      cpat.out,
      tcode.out,
      stopfree.out,
    ) \
    | set { analyses_done }

    qc_analyze(start, analyses_done)
}

workflow {
  main:
    analyze(channel.of('ready'))

  onComplete:
    try {
      def msg = workflow.success ? "Analyze workflow completed" : "Analyze workflow completed with errors"
      slack_closure(msg)
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e}"
    }
}
