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
    // one parquet per assembly; none at all for a taxon with no transcripts
    tuple val(taxid), path("transcripts/transcripts__*.pq", optional: true)

  script:
  """
  rnac genes utils fetch --taxid ${taxid} --output_dir transcripts
  """
}


// Placeholder for the GNN inference container (infer_genes.py), run once per
// (taxid, assembly). It must write assignments.pq with columns region_name, gene_id,
// leaving out regions that are not in a gene (single-member clusters are never genes).
// Fails loudly until it is wired up.
process infer_genes {
  tag { "${taxid} ${assembly}" }
  errorStrategy 'terminate'

  input:
    tuple val(taxid), val(assembly), path(transcripts), path(so_model)

  output:
    tuple val(taxid), val(assembly), path(transcripts), path("assignments.pq")

  script:
  """
  echo "GNN inference is not wired up yet" >&2
  exit 1
  """
}


process make_clusters {
  tag { "${taxid} ${assembly}" }
  memory '16 GB'

  input:
    tuple val(taxid), val(assembly), path(transcripts), path(assignments)

  output:
    // named by assembly: the files of all assemblies are staged together for the plan
    tuple val(taxid), path("clusters__${assembly}.pq")

  script:
  """
  rnac genes utils clusters \
    --transcripts ${transcripts} \
    --assignments ${assignments} \
    --output clusters__${assembly}.pq
  """
}


process fetch_state {
  maxForks 15
  tag { taxid }
  memory '16 GB'

  input:
    val(taxid)

  output:
    tuple val(taxid), path("state.pq"), path("meta.json")

  script:
  """
  rnac genes utils fetch-state --taxid ${taxid} --output state.pq --meta meta.json
  """
}


// Pure: no database access. clusters is empty for a taxon that has lost all its
// transcripts, and then --allow_empty retires its genes.
process plan_genes {
  tag { taxid }
  memory { 32.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    tuple val(taxid), path(state), path(meta), path(clusters), val(allow_empty)

  output:
    tuple val(taxid), path("plan")

  script:
  // a single staged file is a Path, which Groovy iterates by segment: wrap it
  def cluster_args = [clusters].flatten().collect { f -> "--clusters ${f}" }.join(' ')
  def empty_arg = allow_empty ? '--allow_empty' : ''
  """
  rnac genes utils plan \
    --state ${state} \
    --meta ${meta} \
    ${cluster_args} \
    ${empty_arg} \
    --release ${params.release} \
    --out_dir plan
  """
}


process compute_metadata {
  tag { taxid }
  memory '16 GB'

  input:
    tuple val(taxid), path(plan)

  output:
    tuple val(taxid), path(plan), path("metadata.pq")

  script:
  """
  rnac genes utils metadata --genes ${plan}/genes.parquet --taxid ${taxid} --output metadata.pq
  """
}


// One transaction per taxon, one taxon at a time.
process apply_genes {
  tag { taxid }
  maxForks 1
  memory '16 GB'

  input:
    tuple val(taxid), path(plan), path(metadata)

  output:
    val(taxid)

  script:
  """
  rnac genes utils apply \
    --plan_dir ${plan} \
    --metadata ${metadata} \
    --taxid ${taxid} \
    --release ${params.release}
  """
}


workflow genes {
  // _flag is an ordering trigger so this stage can be chained after upstream
  // stages in main.nf. It is not otherwise used: the pipeline is driven by
  // params.genes. Standalone runs pass channel.of(true) (see entry below).
  take: _flag

  main:
  channel.of(true) | fetch_so_model | set { so_model }

  channel.fromPath(params.genes.taxa_query) \
  | fetch_taxids \
  | splitCsv \
  | map { row -> row[0] } \
  | set { taxids }

  taxids | fetch_transcripts | set { fetched }
  taxids | fetch_state | set { state }

  // Transcripts present: cluster each assembly separately with the GNN.
  fetched \
  | filter { _taxid, files -> files } \
  | flatMap { taxid, files ->
      [files].flatten().collect { f ->
        [taxid, f.name.replaceFirst(/^transcripts__/, '').replaceFirst(/\.pq$/, ''), f]
      }
    } \
  | combine( so_model ) \
  | infer_genes \
  | make_clusters \
  | groupTuple \
  | map { taxid, clusters -> [taxid, clusters, false] } \
  | set { with_genes }

  // No transcripts: the plan is empty on purpose and retires the taxon's genes.
  fetched \
  | filter { _taxid, files -> !files } \
  | map { taxid, _files -> [taxid, [], true] } \
  | set { without_genes }

  state \
  | join( with_genes.mix( without_genes ) ) \
  | plan_genes \
  | compute_metadata \
  | apply_genes \
  | set { done }

  // QC: summarise genes created/updated this run once every taxon is stored.
  done | qc_genes

  emit: done
}


workflow {
  genes(channel.of(true))
}
