#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fetch_gene_data {

  input:
    path(query)

  output:
    path("*.pq")

  """
  psql -f $query $PGDATABASE > gene_data.csv
  rnac rf-genes convert gene_data.csv gene_data.pq
  rm gene_data.csv
  """
}

process extract_genes_for_taxid {

  input:
    tuple path(gene_data), val(taxid)

  output:
    tuple path("candidates.pq"), path("transcripts.pq")

  """
  rnac rf-genes extract $gene_data $taxid  --output_genes candidates.pq --output_transcripts transcripts.pq --ensembl_only --pre_gene_set ensembl_genes.pq
  """
}

process preprocess_candidates {

  input:
    tuple path(candidates), path(transcripts)

  output:
    path("features.pq")

  """
  rnac rf-genes preprocess $candidates $transcripts features.pq
  """

}

process split_dataset {

  input:
    path(features)

  output:
    tuple path("train.pq"), path("val.pq"), path("test.pq")

  """
  rnac rf-genes split $features train.pq val.pq test.pq
  """
}


process train_model_nokf{

  input:
    tuple path(train_data), path(val_data), path(test_data)

  output:
    tuple path("transcript_classifier.skops"), path("test_results.md")

  """
  rnac rf-genes train_nokf $train_data transcript_classifier.skops -e comparison
  rnac rf-genes test transcript_classifier.skops $test_data
  """
}


workflow {

  main:
    channel.of("files/random_forest_genes/query.sql") | set { transcript_data_query }

    channel.fromList([10116]) | set { taxids }

    transcript_data_query | fetch_gene_data | set { gene_data }

    gene_data | mix(taxids) | view
    | extract_genes_for_taxid
    | preprocess_candidates
    | split_dataset
    | train_model_nokf
    | set { result }

  publish:
    result >> "results"

}
