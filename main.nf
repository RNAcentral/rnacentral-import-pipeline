#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parse the external database connections file. Previously done inline in
// main.config, but Nextflow 26+ forbids variable declarations in config files.
params.connections = new groovy.json.JsonSlurper().parse(new File(params.connection_file))

// Inferred from params.databases; Nextflow 26 forbids the loop this replaced.
// Expressions live in lib/Utils.groovy so import-data.nf can derive the same.
params.databases.ensembl._any.run = Utils.ensembl_runs(params.databases)
params.needs_publications = params.needs_publications || Utils.needs_publications(params.databases, params.get('skip_publications', false))
params.should_release     = params.should_release || Utils.should_release(params.databases)
params.needs_taxonomy     = params.needs_taxonomy || Utils.needs_taxonomy(params.databases)

include { genes } from './genes'
include { import_data } from './import-data'
include { precompute } from './precompute'
include { analyze } from './analyze'
include { export } from './export'

workflow {
  channel.of('ready') \
  | import_data \
  | ifEmpty('no import') \
  | analyze \
  | ifEmpty('no analysis') \
  | precompute \
  | ifEmpty('no precompute') \
  | genes \
  | ifEmpty('no genes') \
  | export
}
