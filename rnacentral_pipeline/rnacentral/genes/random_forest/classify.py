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
import os
from pathlib import Path

import networkx as nx
import onnxruntime as ort
import polars as pl
from community import community_louvain
from tqdm import tqdm

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

    return_genes = []

    ## try moving the split closer to where it is needed, so we don't store so many tuples
    classifications = classifications.with_columns(
        node_tuple=pl.col("comparison").str.split(" vs ")
    )
    for chromosome_transcripts in tqdm(
        transcripts.partition_by("chromosome"),
        total=transcripts.select("chromosome").unique().height,
    ):
        G = nx.Graph()
        chromosome = chromosome_transcripts.get_column("chromosome").to_list()[0]
        ## partition the classifications per chromosome to make this manageable
        chromosome_classifications = classifications.filter(
            pl.col("chromosome") == chromosome
        )

        # Add all nodes
        node_assembly_lookup = {
            n: a
            for n, a in zip(
                chromosome_transcripts.get_column("region_name").unique().to_list(),
                chromosome_transcripts.get_column("assembly_id").to_list(),
            )
        }
        nodes = [
            (n, {"assembly": a})
            for n, a in zip(
                chromosome_transcripts.get_column("region_name").to_list(),
                chromosome_transcripts.get_column("assembly_id").to_list(),
            )
        ]

        G.add_nodes_from(nodes)

        # Add weighted edges
        edges = (
            chromosome_classifications.filter(pl.col("prediction") == 1)
            .select(["node_tuple", "probability"])
            .collect()
        )

        edge_tuples = []
        for idx, row in enumerate(edges.iter_rows()):
            source, target = row[0][0], row[0][1]
            probability = row[1]
            edge_tuples.append((source, target, probability))
            if idx > 0 and idx % 100_000 == 0:
                G.add_weighted_edges_from(edge_tuples)
                edge_tuples = []
        if len(edge_tuples) > 0:
            G.add_weighted_edges_from(edge_tuples)
        print("edges added")
        # Find communities (genes) using Louvain method
        partition = community_louvain.best_partition(G)

        # Group transcripts by community
        communities = {}
        for transcript, community_id in partition.items():
            if community_id not in communities:
                communities[community_id] = set()
            communities[community_id].add(transcript)
        genes = list(communities.values())
        for g in genes:
            if node_assembly_lookup.get(list(g)[0], None) is not None:
                return_genes.append((g, node_assembly_lookup[list(g)[0]]))

    return return_genes


def run_classification(model_path, features):
    """
    Run the model over the features and classify every single pair
    """
    excluded_columns = ["comparison", "label"]

    sess_opt = ort.SessionOptions()
    sess_opt.execution_mode = ort.ExecutionMode.ORT_PARALLEL
    if os.getenv("SLURM_JOB_ID") is not None:
        cpus_on_node = os.getenv("SLURM_CPUS_ON_NODE") or 1
        sess_opt.intra_op_num_threads = int(cpus_on_node)
        sess_opt.inter_op_num_threads = int(cpus_on_node)

    sess = ort.InferenceSession(
        model_path, providers=["CPUExecutionProvider"], sess_options=sess_opt
    )
    input_name = sess.get_inputs()[0].name  ## gets the probability dict
    label_name = sess.get_outputs()[0].name
    probability_name = sess.get_outputs()[1].name

    classifications = pl.DataFrame(
        {
            "chromosome": [],
            "comparison": [],
            "prediction": [],
            "probability": [],
        },
        schema={
            "chromosome": pl.Utf8,
            "comparison": pl.Utf8,
            "prediction": pl.Int64,
            "probability": pl.Float64,
        },
    )
    if features.height > 10_000:
        chunk_size = 10_000
        for batch_idx in tqdm(
            range(0, features.height, chunk_size),
            total=(features.height // chunk_size) + 1,
        ):
            features_slice = features.slice(batch_idx, chunk_size)

            X = features_slice.select(pl.exclude(excluded_columns)).to_numpy()
            comparisons = features_slice.get_column("comparison").to_numpy()
            features_slice = features_slice.with_columns(
                chromosome=pl.col("comparison")
                .str.split("@")
                .list.get(1)
                .str.split("/")
                .list.first()
            )
            chromosomes = features_slice.get_column("chromosome").to_numpy()
            pred_part, prob_part = sess.run(
                [label_name, probability_name],
                {input_name: X},
            )

            new_classifications = pl.DataFrame(
                {
                    "chromosome": chromosomes,
                    "comparison": comparisons,
                    "prediction": pred_part,
                    "probability": [pr[1] for pr in prob_part],
                },
                schema={
                    "chromosome": pl.Utf8,
                    "comparison": pl.Utf8,
                    "prediction": pl.Int64,
                    "probability": pl.Float64,
                },
            )
            classifications = classifications.vstack(new_classifications)
    else:
        X = features.select(pl.exclude(excluded_columns)).to_numpy()
        comparisons = features.get_column("comparison").to_numpy()
        features = features.with_columns(
            chromosome=pl.col("comparison")
            .str.split("@")
            .list.get(1)
            .str.split("/")
            .list.first()
        )
        chromosomes = features.get_column("chromosome").to_numpy()
        predictions, prob_dict = sess.run(
            [label_name, probability_name], {input_name: X}
        )
        classifications = pl.DataFrame(
            {
                "chromosome": chromosomes,
                "comparison": comparisons,
                "prediction": predictions,
                "probability": [pr[1] for pr in prob_dict],
            },
            schema={
                "chromosome": pl.Utf8,
                "comparison": pl.Utf8,
                "prediction": pl.Int64,
                "probability": pl.Float64,
            },
        )

    ## Rechunk because the many vstacks don't make the memory contiguous, which will slow the next
    ## stage down
    return classifications.rechunk().unique("comparison")


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

    if features.height == 0:
        return pl.DataFrame()

    # Run classification
    classifications = run_classification(model_path, features)

    ## The classifications is the number of edges in the whole genome graph
    ## This is often preposterously huge, meaning the graph construction
    ## is very memory intensive. So, instead we write the classifications
    ## to a parquet file and re-open it in lazy mode. This allows us to do
    ## operations on it without running OOM, at the cost of a little speed.
    classifications.write_parquet("_temp_classifications.parquet")
    del classifications
    classifications = pl.scan_parquet("_temp_classifications.parquet")
    print("Done classifying, now generating genes")
    # Get unique transcript names from features
    transcripts = pl.read_parquet(transcripts_file)

    genes = get_community_genes(classifications, transcripts)
    genes = list(filter(lambda x: len(x[0]) > 1, genes))

    if len(genes) == 0:
        return pl.DataFrame()

    # Name genes
    prefix = data.get_stable_prefix(taxid, conn_str)
    genes_table = data.name_genes(genes, prefix, seed=seed)

    ## Clean up
    os.unlink("_temp_classifications.parquet")

    return genes_table
