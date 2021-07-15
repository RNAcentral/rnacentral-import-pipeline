#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process find_models {
  when { params.cpat.run }

  input:
  val(_flag)

  output:
  path('CPAT-3.0.4/dat/*_logitModel.RData'), emit: rdata
  path('CPAT-3.0.4/dat/*_Hexamer.tsv'), emit: hexamers
  path('cutoffs.csv'), emit: cutoffs

  """
  wget -O cpat.tar.gz 'https://files.pythonhosted.org/packages/c1/3d/de83074cb1b88214db4b48cd894a3ab395d6c5cd1c2c94f14d1950408f6d/CPAT-3.0.4.tar.gz'
  tar xf cpat.tar.gz
  rnac cpat generate-cutoffs CPAT-3.0.4 cutoffs.csv
  """
}

process find_sequences {
  tag { "$model_name" }

  input:
  tuple val(model_name), path(rdata), path(hexamer), val(taxid), path(query)

  output:
  tuple val(model_name), path(rdata), path(hexamer), path('sequences/*.fasta')

  """
  psql -v ON_ERROR_STOP=1 -v "taxid=$taxid" -f "$query" "$PGDATABASE" > raw.json
  mkdir sequences
  split --lines=${params.cpat.chunk_size} --additional-suffix='.fasta' --filter 'json2fasta - - >> \$FILE' raw.json sequences/$model_name-
  """
}

process cpat_scan {
  tag { "$sequences" }
  container 'rnacentral/cpat'

  input:
  tuple val(model_name), path(data), path(hexamer), path(sequences)

  output:
  tuple val(model_name), path('output.ORF_prob.best.tsv')

  """
  cpat.py -g "$sequences" -d "$data" -x "$hexamer" -o output
  """
}

process parse_results {
  input:
  tuple val(model_name), path('scan-results.tsv'), path(cutoff_info)

  output:
  path('results.csv'), emit: results
  path('orfs.csv'), emit: orfs

  """
  rnac cpat parse $cutoff_info $model_name scan-results.tsv .
  """
}

process store_results {
  input:
  path('results*.csv')
  path('orfs.*.csv')
  path(result_ctl)
  path(orf_ctl)

  """
  split-and-load $result_ctl 'results*.csv' ${params.import_data.chunk_size} cpat-results
  split-and-load $orf_ctl 'orfs*.csv' ${params.import_data.chunk_size} cpat-orfs
  """
}

workflow cpat {
  take: flag
  main:
    Channel.fromPath('files/cpat/results.ctl') | set { load_ctl }
    Channel.fromPath('files/cpat/orfs.ctl') | set { orf_ctl }
    Channel.fromPath('files/cpat/query.sql') | set { query }

    flag | find_models

    find_models.out.rdata \
    | flatten \
    | map { it -> [it.getSimpleName().split("_")[0].toLowerCase(), it] } \
    | set { rdata }

    find_models.out.hexamers \
    | flatten \
    | map { it -> [it.getSimpleName().split("_")[0].toLowerCase(), it] } \
    | set { hexamers }

    rdata \
    | join(hexamers) \
    | map { model_name, rdata, hexamer -> [model_name, rdata, hexamer, params.cpat.taxid_mapping[model_name]] } \
    | combine(query) \
    | find_sequences \
    | flatMap { model_name, rdata, hexamer, seqs -> (seqs instanceof ArrayList) ? seqs.collect { [model_name, rdata, hexamer, it] } : [[model_name, rdata, hexamer, seqs]] } \
    | cpat_scan \
    | combine(find_models.out.cutoffs) \
    | parse_results

    parse_results.out.results | collect | set { data }
    parse_results.out.orfs | collect | set { orfs }

    store_results(data, orfs, load_ctl, orf_ctl)
}

workflow {
  cpat(Channel.of('ready'))
}
