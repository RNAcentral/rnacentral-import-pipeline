# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from pathlib import Path

import networkx as nx
import numpy as np
import onnxruntime as ort
import polars as pl
from community import community_louvain

from rnacentral_pipeline.rnacentral.genes.random_forest import data


def convert_for_load(gene_table, transcripts):
    """
    Convert into the order needed for the existing load ctls

    assembly_id,
    locus_name,
    locus_public_name,
    chromosome,
    strand,
    locus_start,
    locus_stop,
    member_count,
    urs_taxid,
    region_id,
    member_type

    For status:

    assembly_id,
    region_id,
    urs_taxid,
    status

    """
    gene_table = gene_table.with_columns(
        pl.col("members").list.lengths().alias("member_count")
    )
    gene_table = gene_table.with_columns(member_type=pl.lit("member"))
    gene_table = gene_table.with_columns(
        pl.when(pl.col("strand") == "-").then(-1).otherwise(1).alias("strand")
    )
    gene_table = gene_table.explode("members")
    gene_table = gene_table.with_columns(
        urs_taxid=pl.col("members").str.split("@").list.first()
    )

    loadable_table = gene_table.join(
        transcripts.select(["region_name", "region_id"]),
        left_on="members",
        right_on="region_name",
        how="inner",
    )

    loadable_table = loadable_table.select(
        [
            "assembly_id",
            "internal_name",
            "name",
            "chromosome",
            "strand",
            "start",
            "stop",
            "member_count",
            "urs_taxid",
            "region_id",
            "member_type",
        ]
    )

    status_table = loadable_table.select(
        ["assembly_id", "region_id", "urs_taxid"]
    ).with_columns(status=pl.lit("inferred"))

    return loadable_table, status_table


def get_community_genes(classifications, transcripts):
    G = nx.Graph()
    # Add all nodes
    node_assembly_lookup = {
        n: a
        for n, a in zip(
            transcripts.get_column("region_name").unique().to_list(),
            transcripts.get_column("assembly_id").to_list(),
        )
    }
    nodes = [
        (n, {"assembly": a})
        for n, a in zip(
            transcripts.get_column("region_name").unique().to_list(),
            transcripts.get_column("assembly_id").unique().to_list(),
        )
    ]
    G.add_nodes_from(nodes)

    # Add weighted edges
    classifications = classifications.with_columns(
        node_tuple=pl.col("comparison").str.split(" vs ")
    )
    edges = classifications.filter(pl.col("prediction") == 1).select(
        ["node_tuple", "probability"]
    )

    for row in edges.iter_rows():
        source, target = row[0][0], row[0][1]
        probability = row[1]
        G.add_edge(source, target, weight=probability)

    # Find communities (genes) using Louvain method
    partition = community_louvain.best_partition(G)

    # Group transcripts by community
    communities = {}
    for transcript, community_id in partition.items():
        if community_id not in communities:
            communities[community_id] = set()
        communities[community_id].add(transcript)

    genes = list(communities.values())
    genes = [(g, node_assembly_lookup[list(g)[0]]) for g in genes]

    return genes


def run_classification(model_path, features):
    """
    Run the model over the features and classify every single pair
    """
    excluded_columns = ["comparison", "label"]
    X = features.select(pl.exclude(excluded_columns)).to_numpy()

    sess_opt = ort.SessionOptions()
    sess_opt.execution_mode = ort.ExecutionMode.ORT_PARALLEL
    sess_opt.inter_op_num_threads = 1

    sess = ort.InferenceSession(model_path, providers=["CPUExecutionProvider"])
    input_name = sess.get_inputs()[0].name  ## gets the probability dict
    label_name = sess.get_outputs()[0].name
    probability_name = sess.get_outputs()[1].name

    if features.height > 10_000:
        predictions = []
        prob_dict = []
        chunk_size = 1000
        for batch_idx in range(0, features.height, chunk_size):
            pred_part, prob_part = sess.run(
                [label_name, probability_name],
                {input_name: X[batch_idx : batch_idx + chunk_size]},
            )
            predictions.extend(pred_part)
            prob_dict.extend(prob_part)
    else:
        predictions, prob_dict = sess.run(
            [label_name, probability_name], {input_name: X}
        )

    probabilities = [pr[c] for pr, c in zip(prob_dict, predictions)]

    comparisons = features.get_column("comparison").to_numpy()

    classifications = pl.DataFrame(
        {
            "comparison": comparisons,
            "prediction": predictions,
            "probability": probabilities,
        }
    )

    return classifications.unique("comparison")


def run_final_classification(
    features_file, transcripts_file, model_path, taxid, conn_str, seed=20240520
):
    """
    Run the final classification step: load features, apply model, generate genes, and save results.

    Args:
        features_file: Path to the features parquet file
        model_path: Path to the trained model pickle file
        taxid: Taxonomy ID for gene naming
        conn_str: Database connection string
        output_dir: Directory to save the output CSV files
        seed: Random seed for gene naming
    """
    # Validate input files exist
    if not Path(features_file).exists():
        raise FileNotFoundError(f"Features file not found: {features_file}")
    if not Path(model_path).exists():
        raise FileNotFoundError(f"Model file not found: {model_path}")
    if not Path(transcripts_file).exists():
        raise FileNotFoundError(f"Model file not found: {transcripts_file}")

    # Load features
    features = pl.read_parquet(features_file)

    transcripts = pl.read_parquet(transcripts_file)

    # Run classification
    classifications = run_classification(model_path, features)

    # Get unique transcript names from features
    all_comparisons = features.get_column("comparison").to_list()
    transcript_names = set()
    for comparison in all_comparisons:
        parts = comparison.split(" vs ")
        transcript_names.update(parts)

    genes = get_community_genes(classifications, transcripts)
    genes = list(filter(lambda x: len(x[0]) > 1, genes))

    if len(genes) == 0:
        return pl.DataFrame()

    # Name genes
    prefix = data.get_stable_prefix(taxid, conn_str)
    genes_table = data.name_genes(genes, prefix, seed=seed)

    return genes_table
