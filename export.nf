#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { text_search } from './workflows/export/text-search.nf'
include { ftp } from './workflows/export/ftp.nf'
include { sequence_search } from './workflows/export/sequence-search'

include { qc_export } from './workflows/utils/qc'

workflow export {
  take: ready
  main:
    ready | (text_search & ftp & sequence_search) | mix | collect | set { done }

    // Final QC: validate published FTP files + search index comparison.
    done | qc_export
  emit: done
}

workflow {
  export(channel.of('ready'))
}
