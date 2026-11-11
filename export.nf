#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { text_search } from './workflows/export/text-search.nf'
include { ftp } from './workflows/export/ftp.nf'
include { sequence_search } from './workflows/export/sequence-search'
include { sankey } from './workflows/export/sankey'
include { huggingface } from './workflows/export/huggingface'
include { active_sequences } from './workflows/export/active-sequences'
include { database_sequences } from './workflows/export/database-sequences'

include { qc_export } from './workflows/utils/qc'

workflow export {
  take: ready
  main:
    done = channel.empty()

    // One dump of the active sequences feeds both the FTP FASTA and HuggingFace;
    // skip the (very slow) query when neither runs.
    active_json = channel.empty()
    if (params.export.ftp.run || params.export.huggingface.run) {
      active_json = active_sequences()
    }

    // The per-database dumps likewise serve both the FTP and sequence-search exports.
    database_dumps = channel.empty()
    if (params.export.ftp.run || params.export.sequence_search.run) {
      database_dumps = database_sequences()
    }

    if (params.export.search.run) {
      text_search(ready)
      done = done.mix(text_search.out)
    }
    if (params.export.ftp.run) {
      ftp(ready, active_json, database_dumps)
      done = done.mix(ftp.out)
    }
    if (params.export.sequence_search.run) {
      sequence_search(ready, database_dumps)
      done = done.mix(sequence_search.out)
    }
    if (params.export.huggingface.run) {
      huggingface(ready, active_json)
      done = done.mix(huggingface.out)
    }
    if (params.export.sankey.run) {
      sankey(ready)
      done = done.mix(sankey.out)
    }

    // Final QC: validate published FTP files + search index comparison.
    done | collect | qc_export
  emit: done
}

workflow {
  export(channel.of('ready'))
}
