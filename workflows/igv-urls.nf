#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process setup {
  when { params.genome_mapping.run }

  input:
  val(_flag)
  path(query)

  output:
  path('species.csv')

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > species.csv
  """
}

process create_json {
  tag { species }
  errorStrategy 'ignore'

  input:
  tuple val(species), val(assembly), val(taxid), val(division)

  output:
  path("${species}.${assembly}.json")

  """
  set -o pipefail
  rnac genome-mapping igv $species $assembly ${species}.${assembly}.json
  """
}

process merge_json {
  publishDir "${params.export.ftp.publish}/.genome-browser", mode: 'copy'

  input:
  path(json)

  output:
  path('igv_files.json')

  """
  jq -s 'flatten' $json > igv_files.json
  """
}

workflow genome_mapping {
  take: ready
  main:
    Channel.fromPath('files/genome-mapping/find_species.sql').set { find_species }

    setup(ready, find_species) \
    | splitCsv \
    | filter { s, a, t, d -> !params.genome_mapping.species_excluded_from_mapping.contains(s) } \
    | set { genome_info }

    genome_info | create_json | collect | merge_json
}

workflow {
  genome_mapping(Channel.from('ready'))
}
