#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parse the external database connections file. Previously done inline in
// main.config, but Nextflow 26+ forbids variable declarations in config files.
params.connections = new groovy.json.JsonSlurper().parse(new File(params.connection_file))

// Infer needs_publications, should_release, and needs_taxonomy from which
// databases are configured to run. Previously a for loop in nextflow.config,
// which Nextflow 26+ no longer supports in config files.
for (item in params.databases) {
  def db = item.value;
  def will_run = db.get('run', false);
  if (item.key == 'ensembl') {
    def ensembl_will_run = ['vertebrates', 'plants', 'fungi', 'protists', 'metazoa'].any {
      db[it]?.get('run', false)
    }
    will_run = ensembl_will_run
    if (ensembl_will_run) {
      params.databases.ensembl._any.run = ensembl_will_run
    }
  }
  if (will_run) {
    params.needs_publications = true;
    if (db.get('release', true)) {
      params.should_release = true;
    }
    if (db.get('needs_taxonomy', false)) {
      params.needs_taxonomy = true
    }
  }
}

include { genes } from './genes'
include { import_data } from './import-data'
include { precompute } from './precompute'
include { analyze } from './analyze'
include { export } from './export'

workflow {
  Channel.of('ready') \
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
