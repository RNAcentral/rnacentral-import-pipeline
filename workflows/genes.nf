#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fetch_taxids {

  when { params.genes.run }

  input:
    path(taxa_query)

  output:
    path('taxids.csv')

  """
  psql \
    -v ON_ERROR_STOP=1 \
    -f "$taxa_query" \
    "$PGDATABASE" > taxids.csv
  """
}

process fetch_rf_model {
  when { params.genes.run }

  input:
    val(_flag)

  output:
    path("genes_model.pkl")

  """
  wget -O genes_model.pkl ${params.genes.rf_model_url}
  """
}

process fetch_so_model{
  when { params.genes.run }

  input:
    val(_flag)

  output:
    path("so_embedding_model.emb")

  """
  wget -O so_embedding_model.emb ${params.genes.so_model_url}
  """
}



process fetch_transcripts {

  tag { taxid }

  input:
    val(taxid)

  output:
    tuple val(taxid), path("transcripts.pq")


  """
  rnac genes fetch --taxid ${taxid} --output transcripts.pq
  """
}


process preprocess_transcripts {

  tag { taxid }

  input:
    tuple val(taxid), path(transcripts), path(so_model)

  output:
    tuple val(taxid), path('features.pq'), path(transcripts)

  """
  rnac genes preprocess --transcripts_file ${transcripts} --output features.pq --so_model_path ${so_model}
  """
}


process classify_transcripts {

  tag { taxid }

  input:
    tuple val(taxid), path(features), path(transcripts), path(model)

  output:
    tuple val(taxid), path("*.csv")

  """
  rnac genes classify --features_file ${features} --transcripts_file ${transcripts} --model_path ${model} --taxid ${taxid} --output_dir .
  """

}



workflow rnacentral_genes {
  emit: genes

  main:
  Channel.of(true) | fetch_so_model | set { so_model }
  Channel.of(true) | fetch_rf_model | set { rf_model }
  Channel.fromPath(params.genes.taxa_query) | fetch_taxids | set { taxa }

  taxa \
  | splitCsv { row -> row[0] } \
  | fetch_transcripts \
  | combine( so_model ) \
  | preprocess_transcripts \
  | combine( rf_model ) \
  | classify_transcripts \
  | set { genes }

}


workflow {
  rnacentral_genes()
}
