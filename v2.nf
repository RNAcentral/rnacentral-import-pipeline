#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ftp_export } from './ftp-export'
include { genes } from './genes'
include { genome_mapping } from 'genome-mapping'
include { import_data } from './import-data'
include { precompute } from './precompute'
include { qa } from './qa'
include { r2dt } from './r2dt'
include { search_export } from './search-export'
include { sequence_search_export } from './sequence-search-export'

workflow {
  import_data \
  | ifEmpty('no import') \
  | set { imported }

  imported \
  | genome_mapping \
  | ifEmpty('no mapping') \
  | set { genome_mapping_results }

  imported \
  | qa \
  | ifEmpty('no qa') \
  | set { qa_results }

  imported \
  | r2dt \
  | ifEmpty('no r2dt') \
  | set { r2dt_results }

  Channel.empty() \
  | mix(qa_results, genome_mapping_results, r2dt_result) \
  | collect \
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
