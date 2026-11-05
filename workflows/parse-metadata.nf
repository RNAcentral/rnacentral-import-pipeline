include { rfam } from './metadata/rfam'
include { ensembl } from './metadata/ensembl'
include { taxonomy } from './metadata/taxonomy'

workflow parse_metadata {
  main:
    // Taxonomy is always refreshed: any database can carry a taxid NCBI added since
    // the last run, and populate_precompute fails on one missing from rnc_taxonomy.
    channel.empty() \
    | mix(
      rfam(),
      (params.databases.ensembl?.run || params.databases.tarbase?.run || params.databases.lncbase?.run) ? ensembl() : channel.empty(),
      taxonomy(),
    ) \
    | flatten \
    | set { data }
  emit: data
}
