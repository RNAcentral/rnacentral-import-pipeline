#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { cpat } from './workflows/cpat'
include { genome_mapping } from './workflows/genome-mapping'
include { qa } from './workflows/qa'
include { r2dt } from './workflows/r2dt'
include { rfam_scan } from './workflows/rfam-scan'

workflow analyze {
  take: ready
  emit: done
  main:
    ready | (genome_mapping & rfam_scan & r2dt & cpat) | mix | collect | set { done }
}

workflow {
  analyze(Channel.of('ready'))
}
