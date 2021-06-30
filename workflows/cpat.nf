process find_models {
  when { params.cpat.run }

  output:
  path('cpat/dat/*_logitModel.RData'), emit: rdata
  path('cpat/dat/*_Hexamer.tsv'), emit: hexamers
  path('cutoffs.csv'), emit: cutoffs

  """
  wget -O cpat.tar.gz 'https://files.pythonhosted.org/packages/c1/3d/de83074cb1b88214db4b48cd894a3ab395d6c5cd1c2c94f14d1950408f6d/CPAT-3.0.4.tar.gz'
  tar xf cpat.tar.gz/
  rnac cpat generate-cutoffs cpat/ cutoffs.csv
  """
}

process find_sequences {
  tag { "$taxid" }

  input:
  tuple path(rdata), path(hexamer), val(taxid)

  output:
  tuple path(rdata), path(hexamer), path('sequences/*.fasta')

  """
  psql -v ON_ERROR_STOP=1 -v "taxid=$taxid" -f "$query" "$PGDATABASE" > raw.json
  mkdir sequences
  split --lines=${params.cpat.chunk_size} --additional-suffix='.fasta' --filter 'json2fasta - - >> \$FILE' raw.json sequences/$taxid-
  """
}

process cpat_scan {
  tag { "$sequences" }
  container 'rnacentral-import-pipeline/cpat'

  input:
  tuple path(data), path(hexamer), path(sequences)

  output:
  path('output.tsv')

  """
  cpat.py -g "$sequences" -d "$data" -x "$hexamer" -o output.tsv
  """
}

process parse_results {
  input:
  tuple path('scan-results.tsv'), path(cutoff_info)

  output:
  path('results.csv')

  """
  rnac cpat parse $cutoff_info scan-results.tsv results.csv
  """
}

process store_results {
  input:
  path('data*.csv')
  path(load)

  """
  split-and-load $load 'data*.csv' ${params.import_data.chunk_size} cpat-load
  """
}

workflow cpat {
  Channel.fromPath('files/cpat/load.ctl') | set { load_ctl }

  find_models()

  find_models.out.rdata \
  | flatten \
  | map { it -> [it.name.split("_").lower(), it] } \
  | set { rdata }

  find_models.out.hexamers \
  | flatten \
  | map { it -> [it.name.split("_").lower(), it] } \
  | set { hexamers }

  rdata \
  | join(hexamers) |
  | map { model_name, rdata, hexamer -> 
    [rdata, hexamer, params.cdat.taxids[model_name]]
  } \
  | combine(query) \
  | find_sequences \
  | flatMap { rdata, hexamer, seqs ->
      (seqs instanceof ArrayList) ? seqs.collect { [rdata, hexamer, it] } : [[rdata, hexamer, seqs]]
  } \
  | cpat_scan \
  | combine(find_models.out.cutoffs) \
  | parse_results \
  | set { data }

  store_results(data, load_ctl)
}

workflow {
  cpat(Channel.of('ready'))
}
