include { rfam } from './metadata/rfam'
include { ensembl } from './metadata/ensembl'
include { taxonomy } from './metadata/taxonomy'
include { any_ensembl_division_runs } from './utils/ensembl-divisions'

workflow parse_metadata {
  main:
    channel.empty() \
    | mix(
      rfam(),
      any_ensembl_division_runs() ? ensembl() : channel.empty(),
      params.needs_taxonomy ? taxonomy() : channel.empty(),
    ) \
    | flatten \
    | set { data }
  emit: data
}
