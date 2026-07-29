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
  main:
    channel.of("Starting data import pipeline") | slack_message

    channel.empty() \
    | mix(
      parse_databases(),
      parse_metadata(),
    ) \
    | branch { r ->
      terms: r.name == "terms.csv"
      ref_ids: r.name == "ref_ids.csv"
      csv: true
    } \
    | set { results }

    results.terms | batch_lookup_ontology_information | set { term_info }
    results.ref_ids | lookup_ref_ids | set { references }

    results.csv \
    | mix(term_info, references) \
    | load_data \
    | set { post_release }



  emit: post_release
}

workflow {
  main:
    import_data(channel.of('ready'))

  onComplete:
    try {
      slack_closure("Workflow completed ${workflow.success ? 'Ok' : 'with errors'}")
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e}"
    }

  onError:
    try {
      slack_closure("Import pipeline encountered an error and failed")
    } catch (Exception e) {
      log.warn "Could not send Slack notification: ${e}"
    }
}
