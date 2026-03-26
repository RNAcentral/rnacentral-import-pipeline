#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process build_ranges {
  when { params.stopfree.run }

  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.stopfree.db_chunk_size
  """
  rnac upi-ranges --table-name rna $chunk_size ranges.csv
  """
}

process find_sequences {
  when { params.stopfree.run }
  maxForks params.stopfree.query_max_forks

  input:
  tuple val(min), val(max), path(query)

  output:
  path('sequences/*.fasta'), optional: true

  """
  PGOPTIONS='-c max_parallel_workers_per_gather=0' psql -v ON_ERROR_STOP=1 -v "min=$min" -v "max=$max" -f "$query" "$PGDATABASE" > raw.json
  mkdir sequences
  split --lines=${params.stopfree.chunk_size} --additional-suffix='.fasta' --filter '${workflow.launchDir}/bin/json2fasta.py - - >> \$FILE' raw.json sequences/seq-
  """
}

process stopfree_scan {
  tag { "$sequences" }
  maxForks params.stopfree.max_forks

  input:
  path(sequences)

  output:
  path("results.csv")

  """
  PYTHONPATH="${workflow.launchDir}" python "${workflow.launchDir}/rnacentral_pipeline/stopfree/scan.py" "$sequences" . --max-probability ${params.stopfree.max_probability}
  """
}

process store_results {
  when { params.stopfree.load }

  input:
  path('results*.csv')
  path(result_ctl)

  """
  split-and-load $result_ctl 'results*.csv' ${params.import_data.chunk_size} stopfree-results
  """
}

workflow stopfree {
  take: flag
  main:
    if( !params.stopfree.run )
      return

    def query = file(params.stopfree.query)
    def load_ctl = file('files/stopfree/stopfree.ctl')

    def fasta_ch = Channel.of('ready') \
      | build_ranges \
      | splitCsv \
      | map { _table, min, max -> [min, max, query] } \
      | find_sequences \
      | flatMap { seqs -> (seqs instanceof ArrayList) ? seqs : [seqs] }

    fasta_ch \
      | stopfree_scan

    stopfree_scan.out | collect | set { data }

    store_results(data, load_ctl)
}

workflow {
  stopfree(Channel.of('ready'))
}
