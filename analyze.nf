#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { qa } from './workflows/qa'
include { genome_mapping } from './workflows/genome-mapping'
include { r2dt } from './workflows/r2dt'

workflow annotate {
  take: ready
  emit: done
  main:
    ready | (genome_mapping & qa & r2dt) | mix | collect | set { done }
}
