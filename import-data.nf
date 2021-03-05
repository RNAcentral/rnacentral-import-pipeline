#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { lookup_ref_ids } from './workflows/lookup-references'
include { batch_lookup_ontology_information } from './workflows/lookup-ontology-info'
include { parse_databases } from './workflows/parse-databases'
include { parse_metadata } from './workflows/parse-metadata'
include { load_data } from './workflows/load-data'

workflow import_data {
  take: databases
  emit: post_release
  main:
    databases \
    | mix(
      parse_databases(),
      parse_metadata(),
    ) \
    | branch {
      terms: it.name == "terms.csv"
      ref_ids: it.name == "ref_ids.csv"
      csv: true
    } \
    | set { results }

    results.terms | batch_lookup_ontology_information | set { term_info }
    results.ref_ids | lookup_ref_ids | set { references }

    results.csv \
    | mix(term_info, references) \
    | load_data \
    | set { post_release }
}

workflow {
  import_data()
}
