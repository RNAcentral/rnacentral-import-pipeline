#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { search_new_ids } from './workflows/litscan/litscan-search-new-ids'
include { load_organisms } from './workflows/litscan/litscan-load-organisms'
include { find_organisms } from './workflows/litscan/litscan-find-organisms'
include { find_retracted_articles } from './workflows/litscan/litscan-retracted-articles'
include { find_manually_annotated } from './workflows/litscan/litscan-manually-annotated'
include { export_articles } from './workflows/litscan/litscan-export-articles'
include { export_metadata } from './workflows/litscan/litscan-export-metadata'

workflow {
  search_new_ids \
  | load_organisms \
  | find_organisms \
  | find_retracted_articles \
  | find_manually_annotated \
  | export_articles \
  | export_metadata
}
