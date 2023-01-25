#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { search_new_ids } from './references-search-new-ids'
include { load_organisms } from './references-load-organisms'
include { find_organisms } from './references-find-organisms'
include { manually_annotated } from './references-manually-annotated'
include { export_metadata } from './references-export-metadata'
include { export_articles } from './references-export-articles'

workflow {
  search_new_ids \
  | load_organisms \
  | find_organisms \
  | manually_annotated \
  | export_articles \
  | export_metadata
}
