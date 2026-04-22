#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { search_new_ids }       from './workflows/litscan/litscan-search-new-ids'
include { scan_new_ids }         from './workflows/litscan/litscan-scan-new-ids'
include { load_organisms }       from './workflows/litscan/litscan-load-organisms'
include { find_organisms }       from './workflows/litscan/litscan-find-organisms'
include { find_retracted_articles } from './workflows/litscan/litscan-retracted-articles'
include { find_manually_annotated } from './workflows/litscan/litscan-manually-annotated'
include { export_articles }      from './workflows/litscan/litscan-export-articles'
include { export_metadata }      from './workflows/litscan/litscan-export-metadata'

workflow {
  // 1. Discover new IDs and register them in litscan_job (status='pending').
  //    Emits a channel of files listing the newly registered IDs.
  search_new_ids()

  // 2. Scan each new ID against Europe PMC (one SLURM job per ID).
  //    Waits for all scans to complete before proceeding.
  scan_new_ids(search_new_ids.out)

  // 3. Downstream export stages (unchanged from previous pipeline).
  load_organisms(scan_new_ids.out)
  find_organisms(load_organisms.out)
  find_retracted_articles(find_organisms.out)
  find_manually_annotated(find_retracted_articles.out)
  export_articles(find_manually_annotated.out)
  export_metadata(find_manually_annotated.out, export_articles.out)
}
