#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { lookup_ref_ids } from './workflows/lookup-references'
include { batch_lookup_ontology_information } from './workflows/lookup-ontology-info'
include { parse_databases } from './workflows/parse-databases'
include { parse_metadata } from './workflows/parse-metadata'
include { load_data } from './workflows/load-data'
include { qc_import } from './workflows/utils/qc'
include { slack_message } from './workflows/utils/slack'
include { slack_closure } from './workflows/utils/slack'

// Promote a delta parse's manifest into rnc_import_manifest, once the release has
// committed. Gated on should_release so a run that does not release never advances
// the manifest ahead of the database. See docs/incremental-parsing.md.
process apply_manifest {
  cache false
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  when: { params.get('should_release', false) }

  input:
  tuple path(manifest), val(_ready)

  output:
  val('done')

  script:
  """
  rnac manifest apply $manifest
  """
}

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
      terms: r.name == "terms.csv" || r.name == "terms.parquet"
      ref_ids: r.name == "ref_ids.csv" || r.name == "ref_ids.parquet"
      manifest: r.name == "manifest.csv"
      csv: true
    } \
    | set { results }

    results.terms | batch_lookup_ontology_information | set { term_info }
    results.ref_ids | lookup_ref_ids | set { references }

    results.csv \
    | mix(term_info, references) \
    | load_data \
    | set { post_release }

    // Final import step: QC — per-database rows imported this release.
    post_release | qc_import

    // deletions.csv rides the normal csv stream into load_data (staged via
    // deletions.ctl). manifest.csv is applied only after the release completes.
    results.manifest \
    | combine(post_release) \
    | apply_manifest

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
