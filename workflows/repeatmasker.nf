#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process build_ranges {
  when: { params.repeatmasker?.run }

  input:
  val(_flag)

  output:
  path('ranges.csv')

  script:
  def chunk_size = params.repeatmasker.db_chunk_size
  """
  rnac upi-ranges --table-name rna $chunk_size ranges.csv
  """
}

process find_sequences {
  when: { params.repeatmasker?.run }

  input:
  tuple val(min), val(max), path(query)

  output:
  path('sequences/*.fasta'), optional: true

  script:
  """
  psql -v ON_ERROR_STOP=1 -v "min=$min" -v "max=$max" -v "min_len=${params.repeatmasker.min_len}" -f "$query" "\$PGDATABASE" > raw.json
  mkdir sequences
  split --lines=${params.repeatmasker.chunk_size} --additional-suffix='.fasta' --filter '${workflow.launchDir}/bin/json2fasta.py - - >> \$FILE' raw.json sequences/seq-
  """
}

process repeatmasker_scan {
  tag { "$sequences" }
  container "${params.repeatmasker.container}"
  containerOptions params.repeatmasker.container_options

  input:
  path(sequences)

  output:
  tuple path(sequences), path("${sequences}.out")

  script:
  """
  RepeatMasker -engine hmmer -species root -pa ${params.repeatmasker.cpus} -dir . $sequences
  # RepeatMasker omits the .out file when a chunk has zero hits; create an empty
  # one so the parser still emits absence rows for every scanned sequence.
  touch ${sequences}.out
  """
}

process parse_results {
  input:
  tuple path(sequences), path(repeatmasker_out)

  output:
  path('results.csv'), emit: results
  path('features.csv'), emit: features

  script:
  """
  rnac repeatmasker parse $sequences $repeatmasker_out .
  """
}

process store_results {
  when: { params.repeatmasker?.load }

  memory 9.GB

  input:
  path('repeatmasker-results*.csv')
  path('repeatmasker-features*.csv')
  path(result_ctl)
  path(feature_ctl)

  script:
  """
  split-and-load $result_ctl 'repeatmasker-results*.csv' ${params.import_data.chunk_size} repeatmasker-results
  split-and-load $feature_ctl 'repeatmasker-features*.csv' ${params.import_data.chunk_size} repeatmasker-features
  """
}

workflow repeatmasker {
  take: flag
  main:
    if( !params.repeatmasker.run ) {
      Channel.of('repeatmasker skipped') | set { done }
    } else {

    def query = file(params.repeatmasker.query)
    def result_ctl = file('files/repeatmasker/repeatmasker.ctl')
    def feature_ctl = file('files/repeatmasker/features.ctl')

    Channel.of('ready') \
      | build_ranges \
      | splitCsv \
      | map { _table, min, max -> [min, max, query] } \
      | find_sequences \
      | flatMap { seqs -> (seqs instanceof ArrayList) ? seqs : [seqs] } \
      | repeatmasker_scan \
      | parse_results

    parse_results.out.results | collect | set { results }
    parse_results.out.features | collect | set { features }

    store_results(results, features, result_ctl, feature_ctl)
    results | map { _ -> 'repeatmasker done' } | first | set { done }
    }
  emit: done
}

workflow {
  repeatmasker(Channel.of('ready'))
}
