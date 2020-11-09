include { five_s_rrnadb } from './databases/5srrnadb'
include { gtrnadb } from './databases/gtrnadb'
include { ensembl } from './databases/ensembl'

workflow parse_databases {
  emit: data
  main:
    Channel.empty()
    | mix(
      five_s_rrnadb(),
      gtrnadb(),
    ) \
    | flatten \
    | set { data }
}
