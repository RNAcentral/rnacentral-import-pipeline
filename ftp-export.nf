#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { export_coordinates } from './workflows/ftp/coordinates'
include { id_mapping } from './workflows/ftp/id-mapping'
include { ensembl_export } from './workflows/ftp/ensembl'
include { fasta_export } from './workflows/ftp/sequences'

process release_note {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  publishDir "${params.ftp_export.publish}/", mode: 'copy'
  when: params.ftp_export.release_note.run

  input:
  path(template_file)

  output:
  path('release_notes.txt')

  """
  rnac ftp-export release-note ${template_file} ${params.release} release_notes.txt
  mkdir -p ${params.ftp_export.publish}/help-requests
  """
}

process md5 {
  publishDir "${params.ftp_export.publish}/md5/", mode: 'copy'
  container ''

  when: params.ftp_export.md5.run

  input:
  path(query)
  path('template.txt')

  output:
  path("example.txt")
  path("md5.tsv.gz")
  path("readme.txt")

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > md5.tsv
  head -10 md5.tsv > example.txt
  gzip md5.tsv
  cat template.txt > readme.txt
  """
}

process rfam_annotations {
  publishDir "${params.ftp_export.publish}/rfam/", mode: 'copy'
  when: params.ftp_export.rfam.annotations.run

  input:
  path(query)
  path('template.txt')

  output:
  path('rfam_annotations.tsv.gz')
  path("example.txt")
  path("readme.txt")

  """
  psql -v ON_ERROR_STOP=1 -f "$query" "$PGDATABASE" > rfam_annotations.tsv
  head rfam_annotations.tsv > example.txt
  gzip rfam_annotations.tsv
  cat template.txt > readme.txt
  """
}

process rfam_go_matches {
  memory params.ftp_export.rfam.go_annotations.memory
  publishDir "${params.ftp_export.publish}/go_annotations/", mode: 'copy'
  when: params.ftp_export.rfam.go_annotations.run

  input:
  path(query)

  output:
  path("rnacentral_rfam_annotations.tsv.gz")

  """
  psql -f "$query" "$PGDATABASE" > raw_go.json
  rnac ftp-export rfam-go-annotations raw_go.json rnacentral_rfam_annotations.tsv
  gzip rnacentral_rfam_annotations.tsv
  """
}

process gpi {
  memory params.ftp_export.gpi.memory
  publishDir "${params.ftp_export.publish}/gpi/", mode: 'copy'

  output:
  path("rnacentral.gpi*")

  """
  rnac ftp-export gpi rnacentral.gpi
  gzip -k rnacentral.gpi
  """
}

workflow ftp_export {

  Channel.fromPath('files/ftp-export/md5/md5.sql') | set { md5_query }
  Channel.fromPath('files/ftp-export/md5/readme.txt') | set { md5_template }
  md5(md5_query, md5_template)

  Channel.fromPath('files/ftp-export/release_note.txt') | release_note
  Channel.fromPath('files/ftp-export/go_annotations/rnacentral_rfam_annotations.sql') | rfam_go_matches

  Channel.fromPath('files/ftp-export/rfam/rfam-annotations.sql') | set { rfam_annotation_query }
  Channel.fromPath('files/ftp-export/rfam/readme.txt') | set { rfam_readme }
  rfam_annotations(rfam_annotation_query, rfam_readme)

  gpi()
  id_mapping()
  export_coordinates()
  ensembl_export()
  fasta_export()
}

workflow {
  ftp_export()
}
