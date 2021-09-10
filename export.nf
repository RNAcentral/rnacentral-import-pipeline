#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { text_search } from './workflows/export/text-search'
include { ftp } from './workflows/export/ftp'
include { sequence_search } from './workflows/export/sequence-search'

workflow export {
  take: ready
  emit: done
  main:
    ready | (text_search & ftp & sequence_search) | mix | collect | set { done }
}

workflow {
  export(Channel.of('ready'))
}
