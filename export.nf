#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { text_search } from './workflows/export/text-search.nf'
include { ftp } from './workflows/export/ftp.nf'
include { active_sequences } from './workflows/export/active-sequences'
include { sequence_search } from './workflows/export/sequence-search'

include { qc_export } from './workflows/utils/qc'

workflow export {
  take: ready
  main:
    text_search(ready)
    ftp(ready, active_sequences())
    sequence_search(ready)
    text_search.out | mix(sequence_search.out) | collect | set { done }

    // Final QC: validate published FTP files + search index comparison.
    done | qc_export
  emit: done
}

workflow {
  export(channel.of('ready'))
}
