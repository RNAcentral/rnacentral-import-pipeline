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
    tuple val(taxid), path("transcripts/transcripts__*.pq", arity: '0..*')

  script:
  """
  rnac genes utils fetch --taxid ${taxid} --output_dir transcripts
  """
}


// GNN inference (ghcr.io/rnacentral/rnacentral-gene-gnn-inference) on a GPU. Each job
// takes a batch of transcripts__<taxid>__<assembly>.pq files (many taxa) and writes
// assignments__<taxid>__<assembly>.pq (region_name, gene_id) for each; regions not in
// a gene are left out, since single-member clusters are never genes.
// Mirrors gnn_genes in rnacentral-genes-history/main_replay.nf, which built the history.
process infer_genes {
  tag { "${transcripts.size()} files" }
  container params.genes.gnn.image
  // the image is rooted at /app and /opt, so bind the launch dir only (cf. r2dt_container)
  containerOptions "--contain --workdir $baseDir/work/tmp --bind $baseDir --nv"
  // the accelerator directive gave no --gres on Codon (jobs landed on CPU nodes), so ask SLURM directly
  clusterOptions "--gres=gpu:${params.genes.gnn.gpu_type}:1"
  cpus params.genes.gnn.cpus       // window feature building is threaded
  memory { "${params.genes.gnn.memory_gb * task.attempt} GB" }
  time { "${params.genes.gnn.hours * task.attempt}h" }  // for the whole batch
  // ponytail: a retry reruns the whole batch; skip finished files if retries get costly
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 2

  input:
    path(transcripts)
    path(so_model)

  output:
    path("assignments__*.pq")

  script:
  """
  set -euo pipefail
  for f in transcripts__*.pq; do
    echo "== \$f"
    GNN_CHECKPOINT=${params.genes.gnn.checkpoint} rnacentral-gene-gnn \
      --input "\$f" \
      --rnac-so-embeddings ${so_model} \
      --device gpu \
      --dtype ${params.genes.gnn.dtype} \
      --threshold ${params.genes.gnn.threshold} \
      --output "assignments__\${f#transcripts__}"
  done
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

  // Transcripts present: one file per (taxid, assembly), keyed by "<taxid>__<assembly>".
  fetched \
  | filter { _taxid, files -> files } \
  | flatMap { taxid, files ->
      [files].flatten().collect { f ->
        def key = f.name.replaceFirst(/^transcripts__/, '').replaceFirst(/\.pq$/, '')
        [key, taxid, key.substring("${taxid}__".length()), f]
      }
    } \
  | set { transcripts }

  // Pool every taxon's files into a fixed number of GPU jobs, largest file first onto
  // the least-loaded batch, so one big genome does not leave the other GPUs idle.
  transcripts \
  | map { _key, _taxid, _assembly, f -> f } \
  | collect \
  | flatMap { files ->
      def n = params.genes.gnn.gpu_jobs as int
      def batches = (0..<n).collect { [] }
      def load = [0L] * n
      files.sort(false) { f -> -f.size() }.each { f ->
        def i = load.indexOf(load.min())
        batches[i] << f
        load[i] += f.size()
      }
      batches.findAll { b -> b }
    } \
  | set { batches }

  infer_genes(batches, so_model.first())

  // Back to one (taxid, assembly) per item; a file without assignments is an error.
  infer_genes.out \
  | flatten \
  | map { a -> [a.name.replaceFirst(/^assignments__/, '').replaceFirst(/\.pq$/, ''), a] } \
  | set { assigned }

  transcripts \
  | join( assigned, failOnMismatch: true ) \
  | map { _key, taxid, assembly, f, a -> [taxid, assembly, f, a] } \
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
