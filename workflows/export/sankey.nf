#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process sankey {
  publishDir "${params.export.sankey.publish}/", mode: 'copy'

  input:
  val(_ready)

  output:
  path('*.png')

  when: params.export.sankey.run

  script:
  """
  rnac sankey .
  """
}
