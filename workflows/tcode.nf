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
  split --lines=${params.tcode.chunk_size} --additional-suffix='.fasta' --filter 'json2fasta - - >> \$FILE' raw.json sequences/seq-
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
  # Extract the first FASTA ID for fallback output.
  seq_id=\$(awk '/^>/{print substr(\$0,2); exit}' ${sequences} | awk '{print \$1}')
  seq_len=0

  # If TCODE fails, write a minimal output so parsing yields null scores.
  if ! tcode -sequence ${sequences} -outfile ${sequences.simpleName}.tcode.out -window ${params.tcode.window}; then
    printf "# Sequence: %s\\n# Total_length: %s\\n" "\$seq_id" "\$seq_len" > ${sequences.simpleName}.tcode.out
    exit 0
  fi
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
