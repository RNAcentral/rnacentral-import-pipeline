#!/usr/bin/env nextflow

include { rfam_scan } from './workflows/rfam-scan'

workflow qa {
  emit: flag
  main:
    rfam() | set { flag }
}

workflow {
  qa()
}
