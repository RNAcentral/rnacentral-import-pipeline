#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { search_new_ids } from './litscan-search-new-ids'
include { load_organisms } from './litscan-load-organisms'
include { find_organisms } from './litscan-find-organisms'
include { find_retracted_articles } from './litscan-retracted-article'
include { find_manually_annotated } from './litscan-manually-annotated'
include { export_articles } from './litscan-export-articles'
include { export_metadata } from './litscan-export-metadata'

workflow {
  search_new_ids \
  | load_organisms \
  | find_organisms \
  | find_retracted_articles \
  | find_manually_annotated \
  | export_articles \
  | export_metadata
}
