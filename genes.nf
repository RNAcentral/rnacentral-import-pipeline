#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { qc_genes } from './workflows/utils/qc'

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

process fetch_so_model{
  input:
    val(_flag)

  output:
    path("so_embedding_model.parquet")

  script:
  """
  wget -O so_embedding_model.parquet ${params.genes.so_model_url}
  """
}



process fetch_transcripts {
  memory '32 GB'
  maxForks 15
  tag { taxid }

  input:
    val(taxid)

  output:
    tuple val(taxid), path("transcripts.pq")


  script:
  """
  rnac genes utils fetch --taxid ${taxid} --output transcripts.pq
  """
}


// Placeholder for the GNN inference container (infer_genes.py). It must emit
// genes_<taxid>.json in the name_genes shape; fails loudly until it is wired up.
process infer_genes {
  tag { taxid }

  input:
    tuple val(taxid), path(transcripts), path(so_model)

  output:
    tuple val(taxid), path("genes_${taxid}.json")

  script:
  """
  echo "GNN inference is not wired up yet" >&2
  exit 1
  """
}


process fetch_previous_genes {
  maxForks 15
  tag { taxid }
  memory '8 GB'

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
  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

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
  memory 4.GB
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
  memory '16 GB'
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
  // params.genes. Standalone runs pass channel.of(true) (see entry below).
  take: _flag

  main:
  channel.of(true) | fetch_so_model | set { so_model }
  channel.fromPath(params.genes.taxa_query) | fetch_taxids | set { taxa }

  taxa \
  | splitCsv \
  | map { row -> row[0] } \
  | fetch_transcripts \
  | filter {_taxid, transcripts -> transcripts.size() > 0 } \
  | combine( so_model ) \
  | infer_genes \
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

  // QC: summarise genes created/updated this run once every taxon is stored.
  done | qc_genes

  emit: done
}


workflow {
  genes(channel.of(true))
}
