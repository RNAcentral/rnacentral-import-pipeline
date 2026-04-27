#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { lookup_ref_ids } from './workflows/lookup-references'
include { batch_lookup_ontology_information } from './workflows/lookup-ontology-info'
include { parse_databases } from './workflows/parse-databases'
include { parse_metadata } from './workflows/parse-metadata'
include { load_data } from './workflows/load-data'
include { slack_message } from './workflows/utils/slack'
include { slack_closure } from './workflows/utils/slack'

workflow import_data {
  take: _flag
  emit: post_release
  main:
    Channel.of("Starting data import pipeline") | slack_message

    Channel.empty() \
    | mix(
      parse_databases(),
      parse_metadata(),
    ) \
    | branch {
      // Match either CSV or Parquet emissions of these auxiliary tables —
      // parsers pick one format based on params.writer_format. Both branches
      // route to bespoke processing (ontology lookup / publication lookup);
      // their outputs feed back into the load stream.
      terms: it.name == "terms.csv" || it.name == "terms.parquet"
      ref_ids: it.name == "ref_ids.csv" || it.name == "ref_ids.parquet"
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
  import_data(Channel.of('ready'))
}

workflow.onError {
  slack_closure("Import pipeline encountered an error and failed")
}

workflow.onComplete {
  slack_closure("Workflow completed ${workflow.success ? 'Ok' : 'with errors'}")
}
