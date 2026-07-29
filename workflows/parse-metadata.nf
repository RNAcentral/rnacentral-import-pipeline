include { rfam } from './metadata/rfam'
include { ensembl } from './metadata/ensembl'
include { taxonomy } from './metadata/taxonomy'

workflow parse_metadata {
  main:
    channel.empty() \
    | mix(
      rfam(),
      params.databases.ensembl._any.run ? ensembl() : channel.empty(),
      params.needs_taxonomy ? taxonomy() : channel.empty(),
    ) \
    | flatten \
    | set { data }
  emit: data
}
