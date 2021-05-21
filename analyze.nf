#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { qa } from './workflows/qa'
include { genome_mapping } from './workflows/genome-mapping'

workflow annotate {
  take: ready
  emit: done
  main:
    ready | (genome_mapping & qa) | mix | collect | set { done }
}
