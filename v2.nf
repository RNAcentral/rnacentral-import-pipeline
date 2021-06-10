#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ftp_export } from './ftp-export'
include { genes } from './genes'
include { import_data } from './import-data'
include { precompute } from './precompute'
include { analyze } from './analyze'
include { search_export } from './search-export'
include { sequence_search_export } from './sequence-search-export'

workflow {
  import_data \
  | ifEmpty('no import') \
  | set { imported }

  imported \
  | analyze \
  | ifEmpty('no analysis') \
  | set { analyze_results }

  Channel.empty() \
  | mix(analyze_results) \
  | map { 'done' } \
  | precompute \
  | ifEmpty('no precompute') \
  | genes \
  | ifEmpty('no genes') \
  | set { post_precompute }

  post_precompute | search_export
  post_precompute | sequence_search_export
  post_precompute | ftp_export
}
