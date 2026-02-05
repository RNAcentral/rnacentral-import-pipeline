#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process build_ranges {
  when { params.tcode.run }

  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.tcode.db_chunk_size
  """
  rnac upi-ranges --table-name rna $chunk_size ranges.csv
  """
}

process find_sequences {
  when { params.tcode.run }

  input:
  tuple val(min), val(max), path(query)

  output:
  path('sequences/*.fasta'), optional: true

  """
  psql -v ON_ERROR_STOP=1 -v "min=$min" -v "max=$max" -v "min_len=${params.tcode.min_len}" -f "$query" "$PGDATABASE" > raw.json
  mkdir sequences
  split --lines=${params.tcode.chunk_size} --additional-suffix='.fasta' --filter '${workflow.launchDir}/bin/json2fasta.py - - >> \$FILE' raw.json sequences/seq-
  """
}

process tcode_scan {
  tag { "$sequences" }
  container 'biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1'

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
  rnac tcode parse $tcode_out .
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

    def query = file(params.tcode.query)
    def load_ctl = file('files/tcode/tcode.ctl')

    def fasta_ch = Channel.of('ready') \
      | build_ranges \
      | splitCsv \
      | map { _table, min, max -> [min, max, query] } \
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
