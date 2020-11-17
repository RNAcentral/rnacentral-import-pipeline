include { rfam } from './metadata/rfam'
include { ensembl } from './metadata/ensembl'
include { taxonomy } from './metadata/taxonomy'

workflow parse_metadata {
  emit: data
  main:
    Channel.empty() \
    | mix(
      rfam(), 
      ensembl(),
      taxonomy(),
    ) \
    | flatten \
    | set { data }
}
