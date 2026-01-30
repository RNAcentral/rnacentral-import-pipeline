#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_sequences {
  when { params.tcode.run }

  input:
  tuple val(taxid), path(query)

  output:
  path('sequences/*.fasta')

  """
  psql -v ON_ERROR_STOP=1 -v "taxid=$taxid" -f "$query" "$PGDATABASE" > raw.json
  mkdir sequences
  split --lines=${params.tcode.chunk_size} --additional-suffix='.fasta' --filter '${workflow.launchDir}/bin/json2fasta.py - - >> \$FILE' raw.json sequences/seq-
  """
}

process tcode_scan {
  tag { "$sequences" }
  container 'biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1'
  containerOptions '--platform linux/amd64'

  input:
  path(sequences)

  output:
  path("${sequences.simpleName}.tcode.out")

  """
  tcode -sequence ${sequences} -outfile ${sequences.simpleName}.tcode.out -window ${params.tcode.window}
  """
}

process parse_results {
  input:
  path(tcode_out)

  output:
  path("results.csv")

  """
  ${workflow.launchDir}/bin/rnac tcode parse $tcode_out .
  """
}

process store_results {
  when { params.tcode.load }

  input:
  path('results*.csv')
  path(result_ctl)

  """
  split-and-load $result_ctl 'results*.csv' ${params.import_data.chunk_size} tcode-results
  """
}

workflow tcode {
  take: flag
  main:
    if( !params.tcode.run )
      return

    Channel.fromPath(params.tcode.query) | set { query }
    Channel.fromPath('files/tcode/tcode.ctl') | set { load_ctl }

    def fasta_ch = Channel.of(params.tcode.taxid) \
      | combine(query) \
      | find_sequences \
      | flatMap { seqs -> (seqs instanceof ArrayList) ? seqs : [seqs] }

    fasta_ch \
      | tcode_scan \
      | parse_results

    parse_results.out | collect | set { data }

    store_results(data, load_ctl)
}

workflow {
  tcode(Channel.of('ready'))
}
