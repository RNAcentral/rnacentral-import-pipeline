#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fetch_taxids {

  input:
    path(taxa_query)

  output:
    path('taxids.csv')

  script:
  """
  psql \
    -v ON_ERROR_STOP=1 \
    -f "${taxa_query}" \
    "\$PGDATABASE" > taxids.csv
  """
}

process fetch_rf_model {
  input:
    val(_flag)

  output:
    path("genes_model.pkl")

  script:
  """
  wget -O genes_model.pkl ${params.genes.rf_model_url}
  """
}

process fetch_so_model{
  input:
    val(_flag)

  output:
    path("so_embedding_model.emb")

  script:
  """
  wget -O so_embedding_model.emb ${params.genes.so_model_url}
  """
}



process fetch_transcripts {
  maxForks 15
  tag { taxid }

  input:
    val(taxid)

  output:
    tuple val(taxid), path("transcripts.pq")


  script:
  """
  rnac genes infer fetch --taxid ${taxid} --output transcripts.pq
  """
}


process preprocess_transcripts {
  memory '64 GB'
  tag { taxid }

  input:
    tuple val(taxid), path(transcripts), path(so_model)

  output:
    tuple val(taxid), path('features.pq'), path(transcripts)

  script:
  """
  rnac genes infer preprocess --transcripts_file ${transcripts} --output features.pq --so_model_path ${so_model}
  """
}


process classify_transcripts {
  memory '64 GB'
  tag { taxid }

  input:
    tuple val(taxid), path(features), path(transcripts), path(model)

  output:
    tuple val(taxid), path("*.json")

  script:
  """
  rnac genes infer classify --features_file ${features} --transcripts_file ${transcripts} --model_path ${model} --taxid ${taxid} --output_dir .
  """

}


process fetch_previous_genes {
  maxForks 15
  tag { taxid }

  input:
    tuple val(taxid), path(this_release)

  output:
    tuple val(taxid), path(this_release), path('previous.json'), path('prev_release.txt')

  script:
  """
  rnac genes utils fetch-stored --taxid ${taxid} --output previous.json --release_output prev_release.txt
  """
}


process forward_merge {
  memory '16 GB'
  tag { taxid }

  input:
    tuple val(taxid), path(this_release), path(previous), val(prev_release)

  output:
    tuple val(taxid), path("merged_${taxid}.json")

  script:
  """
  rnac genes utils merge \
    --previous_genes ${previous} \
    --next_genes ${this_release} \
    --output merged_${taxid}.json \
    --prev_release_number ${prev_release} \
    --next_release_number ${params.release}
  """
}


process init_genes {

  tag { taxid }

  input:
    tuple val(taxid), path(this_release)

  output:
    tuple val(taxid), path("merged_${taxid}.json")

  script:
  """
  rnac genes utils init \
    --genes ${this_release} \
    --output merged_${taxid}.json \
    --release_number ${params.release}
  """
}


process store_genes {

  tag { taxid }
  maxForks 1

  input:
    tuple val(taxid), path(merged)

  output:
    tuple val(taxid), path(merged)

  script:
  """
  rnac genes utils store-genes ${merged} --taxid ${taxid}
  """
}


process deactivate_discarded {

  tag { taxid }
  maxForks 1

  input:
    tuple val(taxid), path(merged)

  output:
    tuple val(taxid), path(merged)

  script:
  """
  rnac genes utils deactivate-discarded --merged ${merged} --taxid ${taxid}
  """
}


process process_metadata {

  tag { taxid }

  input:
    tuple val(taxid), path(merged)

  output:
    tuple val(taxid), path("metadata_${taxid}.json")

  script:
  """
  rnac genes utils process-metadata ${merged} metadata_${taxid}.json
  """
}


process store_metadata {

  tag { taxid }
  maxForks 1

  input:
    tuple val(taxid), path(metadata)

  output:
    val(taxid)

  script:
  """
  rnac genes utils store-metadata ${metadata}
  """
}


workflow genes {
  // _flag is an ordering trigger so this stage can be chained after upstream
  // stages in main.nf. It is not otherwise used: the pipeline is driven by
  // params.genes. Standalone runs pass Channel.of(true) (see entry below).
  take: _flag

  main:
  Channel.of(true) | fetch_so_model | set { so_model }
  Channel.of(true) | fetch_rf_model | set { rf_model }
  Channel.fromPath(params.genes.taxa_query) | fetch_taxids | set { taxa }

  taxa \
  | splitCsv \
  | map { row -> row[0] } \
  | fetch_transcripts \
  | filter {_taxid, transcripts -> transcripts.size() > 0 } \
  | combine( so_model ) \
  | preprocess_transcripts \
  | combine( rf_model ) \
  | classify_transcripts \
  | filter {_taxid, genes -> genes.size() > 2 } \
  | fetch_previous_genes \
  | branch { row ->
      existing:  row[3].text.trim() != ''
      new_taxon: true
    } \
  | set { previous_genes }

  // Existing taxon: forward-merge this release onto the stored genes.
  previous_genes.existing \
  | map { taxid, this_release, previous, prev_release ->
      [taxid, this_release, previous, prev_release.text.trim()]
    } \
  | forward_merge \
  | set { merged_existing }

  // New taxon: nothing to merge, just initialise version tracking.
  previous_genes.new_taxon \
  | map { taxid, this_release, _previous, _prev_release -> [taxid, this_release] } \
  | init_genes \
  | set { merged_new }

  // Store, then deactivate genes the merge dropped, then metadata.
  merged_existing.mix( merged_new ) \
  | store_genes \
  | deactivate_discarded \
  | process_metadata \
  | store_metadata \
  | set { done }

  emit: done
}


workflow {
  genes(Channel.of(true))
}
