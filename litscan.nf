#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { search_new_ids } from './litscan-search-new-ids'
include { load_organisms } from './litscan-load-organisms'
include { find_organisms } from './litscan-find-organisms'
include { manually_annotated } from './litscan-manually-annotated'
include { export_metadata } from './litscan-export-metadata'
include { export_articles } from './litscan-export-articles'

workflow {
  search_new_ids \
  | load_organisms \
  | find_organisms \
  | manually_annotated \
  | export_articles \
  | export_metadata
}
