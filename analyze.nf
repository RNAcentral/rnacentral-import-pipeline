#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { qa } from './workflows/qa'
include { genome_mapping } from './workflows/genome-mapping'
include { r2dt } from './workflows/r2dt'
include { cpat } from './workflows/cpat'

workflow analyze {
  take: ready
  emit: done
  main:
    ready | (genome_mapping & qa & r2dt & cpat) | mix | collect | set { done }
}

workflow {
  analyze(Channel.of('ready'))
}
