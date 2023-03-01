#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { genes } from './genes'
include { import_data } from './import-data'
include { precompute } from './precompute'
include { analyze } from './analyze'
include { export } from './export'

workflow {
  import_data \
  | ifEmpty('no import') \
  | analyze \
  | ifEmpty('no analysis') \
  | precompute \
  | ifEmpty('no precompute') \
  | genes \
  | ifEmpty('no genes') \
  | export
}
