#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { export_coordinates } from './ftp/coordinates'
include { id_mapping } from './ftp/id-mapping'
include { ensembl_export } from './ftp/ensembl'
include { fasta_export } from './ftp/sequences'
include { rediportal } from './ftp/rediportal.nf'

process release_note {
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir"

  publishDir "${params.export.ftp.publish}/", mode: 'copy'

  input:
  path(template_file)

  output:
  path('release_notes.txt')

  """
  rnac ftp-export release-note ${template_file} ${params.release} release_notes.txt
  mkdir -p ${params.export.ftp.publish}/help-requests
  """
}

process md5 {
  publishDir "${params.export.ftp.publish}/md5/", mode: 'copy'
  container ''

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
  publishDir "${params.export.ftp.publish}/rfam/", mode: 'copy'

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
  memory params.export.ftp.rfam.go_annotations.memory
  publishDir "${params.export.ftp.publish}/go_annotations/", mode: 'copy'

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
  memory params.export.ftp.gpi.memory
  publishDir "${params.export.ftp.publish}/gpi/", mode: 'copy'

  output:
  path("rnacentral.gpi*")

  """
  rnac ftp-export gpi rnacentral.gpi
  gzip -k rnacentral.gpi
  """
}

workflow ftp {
  take: _flag
  main:
    if (params.export.ftp.run) {
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
      rediportal()
    }
}

workflow {
  ftp(Channel.of('ready'))
}
