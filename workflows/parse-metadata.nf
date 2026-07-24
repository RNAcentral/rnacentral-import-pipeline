include { rfam } from './metadata/rfam'
include { ensembl } from './metadata/ensembl'
include { taxonomy } from './metadata/taxonomy'

workflow parse_metadata {
  main:
    Channel.empty() \
    | mix(
      rfam(),
      params.databases.ensembl._any.run ? ensembl() : Channel.empty(),
      params.needs_taxonomy ? taxonomy() : Channel.empty(),
    ) \
    | flatten \
    | set { data }
  emit: data
}
