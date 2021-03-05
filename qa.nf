#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { rfam_scan } from './workflows/rfam-scan'

workflow qa {
  take: ready
  emit: done
  main:
    ready | rfam_scan | set { done }
}

workflow {
  qa(Channel.from('ready'))
}
