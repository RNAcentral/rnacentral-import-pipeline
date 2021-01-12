#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { rfam_scan } from './workflows/rfam-scan'

workflow qa {
  main:
    rfam_scan()
}

workflow for_database {
  take: sequences
  emit: data
    rfam_scan.for_database(sequences)
}

workflow {
  qa()
}
