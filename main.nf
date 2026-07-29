#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parse the external database connections file. Previously done inline in
// main.config, but Nextflow 26+ forbids variable declarations in config files.
params.connections = new groovy.json.JsonSlurper().parse(new File(params.connection_file))

// Whether a database entry is configured to run. Ensembl is special-cased: it
// runs if any of its divisions do.
def willRun(key, db) {
  if (key == 'ensembl') {
    return ['vertebrates', 'plants', 'fungi', 'protists', 'metazoa'].any { d -> db[d]?.get('run', false) }
  }
  return db.get('run', false)
}

// Infer needs_publications, should_release, and needs_taxonomy from which
// databases are configured to run. Previously a for loop in nextflow.config;
// Nextflow 26+ forbids loops in config files and top-level loop statements in
// scripts, so these are derived with collection expressions (param assignments,
// which are permitted at the top level).
params.databases.ensembl._any.run =
  ['vertebrates', 'plants', 'fungi', 'protists', 'metazoa'].any { d -> params.databases.ensembl[d]?.get('run', false) }
params.needs_publications = !params.get('skip_publications', false) && params.databases.any { key, db -> willRun(key, db) }
params.should_release     = params.databases.any { key, db -> willRun(key, db) && db.get('release', true) }
params.needs_taxonomy     = params.databases.any { key, db -> willRun(key, db) && db.get('needs_taxonomy', false) }

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
