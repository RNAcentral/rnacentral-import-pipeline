#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { select } from './workflows/databases/select.nf'
include { import_data } from './import-data'
include { precompute } from './precompute'
include { analyze } from './analyze'
include { export } from './export'

workflow {
  select \
  | import_data \
  | ifEmpty('no import') \
  | analyze \
  | ifEmpty('no analysis') \
  | precompute \
  | ifEmpty('no precompute') \
  | export
}
