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


from functools import lru_cache
from pathlib import Path

## Use rust preprocessing code
## Must be built using the Makefile in an activated environment
import gene_preprocessing as gpp
import numpy as np
import polars as pl
from gensim.models import Word2Vec
from tqdm import tqdm

from rnacentral_pipeline.rnacentral.genes.random_forest import data

empty_features = pl.DataFrame(
    {
        "5p_exon_overlap": [],
        "5p_exon_dta": [],
        "5p_exon_3p_dta": [],
        "exons_overlapping": [],
        "strand": [],
        "type_sim": [],
        "comparison": [],
    },
    schema={
        "5p_exon_overlap": pl.Float64,
        "5p_exon_dta": pl.Int64,
        "5p_exon_3p_dta": pl.Int64,
        "exons_overlapping": pl.Int64,
        "strand": pl.Int64,
        "type_sim": pl.Float32,
        "comparison": pl.Utf8,
    },
)

_SO_MODEL = None


def init_so_model(path: str) -> None:
    global _SO_MODEL
    _SO_MODEL = _load_so_model(path)


def _load_so_model(path: str) -> dict:
    so_model = Word2Vec.load(path)
    so_vec_normalised = {
        key: so_model.wv[key] / np.linalg.norm(so_model.wv[key])
        for key in so_model.wv.key_to_index
    }
    return so_vec_normalised


@lru_cache(maxsize=1024)
def get_type_similarity(type_a: str, type_b: str) -> float:
    if _SO_MODEL is None:
        raise RuntimeError("SO model not initialized - call init_so_model(path) first")
    if type_a > type_b:
        type_a, type_b = type_b, type_a
    return rna_type_similarity(type_a, type_b, _SO_MODEL)


def exon_overlap(exon_a_start, exon_a_end, exon_b_start, exon_b_end):
    """
    Calculates overlap of exon b with exon a
    """
    # Early exit if no overlap possible
    if exon_a_end <= exon_b_start or exon_b_end <= exon_a_start:
        return 0

    length_a = abs(exon_a_end - exon_a_start)
    length_b = abs(exon_b_end - exon_b_start)

    # Early exit for invalid exons
    if length_a == 0 or length_b == 0:
        return 0

    overlap_start = max(exon_a_start, exon_b_start)
    overlap_stop = min(exon_a_end, exon_b_end)
    overlap_length = abs(overlap_stop - overlap_start)

    overlap_a = overlap_length / length_a
    overlap_b = overlap_length / length_b

    return min(overlap_a, overlap_b)


def exon_overlap_tup(exon_a, exon_b):
    # Normalize coordinates
    a_start, a_end = min(exon_a[0], exon_a[1]), max(exon_a[0], exon_a[1])
    b_start, b_end = min(exon_b[0], exon_b[1]), max(exon_b[0], exon_b[1])

    return exon_overlap(a_start, a_end, b_start, b_end)


def distance_2_agreement(exon_a, exon_b):
    """
    Only look at the absolute distance between the start coordinates
    """
    return tuple(np.abs(np.array(exon_a) - np.array(exon_b)))


def rna_type_similarity(so_type_a, so_type_b, so_model):
    """
    Compute dot product similarity between transcript type vectors
    """
    so_vec_a = so_model[so_type_a]
    so_vec_b = so_model[so_type_b]

    sim = np.dot(so_vec_a, so_vec_b)
    return sim


def run_preprocessing(
    transcripts_file,
    regions_data,
    so_model_path,
    nearby_distance,
):
    transcripts = pl.read_parquet(transcripts_file)
    ## Filter out piRNAs
    transcripts = transcripts.filter(pl.col("so_type") != "SO:0001035")
    if transcripts.height > 0:
        # Check if region_ids are present, fetch if needed
        if "region_id" not in transcripts.columns:
            if not Path(regions_data).exists():
                raise ValueError(
                    "Region IDs not found in transcripts file and no database connection provided. "
                    "Please provide --conn_str to fetch region IDs from database."
                )
            transcripts = data.add_assembly_region_ids(
                transcripts, regions_data
            ).select(pl.exclude("assembly_id_right"))

        init_so_model(so_model_path)
        features = polars_preprocessing(transcripts, nearby_distance)
    else:
        features = empty_features.clone()

    return features


def compare_transcripts(transcripts_a, transcripts_b, so_model, label=0):
    """
    Compare two sets of transcripts and generate features. ONly used when creating training data.

    May be removed in future, in favour of a faster preprocessor
    """
    comparisons = set()
    similarity_comparisons = set()
    similarity_cache = {}

    if label is not None:
        features = pl.DataFrame(
            {
                "5p_exon_overlap": [],
                "5p_exon_dta": [],
                "5p_exon_3p_dta": [],
                "exons_overlapping": [],
                "strand": [],
                "type_sim": [],
                "label": [],
                "comparison": [],
            },
            schema={
                "5p_exon_overlap": pl.Float64,
                "5p_exon_dta": pl.Int64,
                "5p_exon_3p_dta": pl.Int64,
                "exons_overlapping": pl.Int64,
                "strand": pl.Int64,
                "type_sim": pl.Float64,
                "label": pl.Int8,
                "comparison": pl.Utf8,
            },
        )
    else:
        features = empty_features.clone()

    for tr_a in transcripts_a.iter_rows(named=True):
        comparison_features = []
        for tr_b in transcripts_b.iter_rows(named=True):
            if tr_a["region_name"] == tr_b["region_name"]:
                continue
            ## Check if we already compared these two transcripts
            normalized_pair = tuple(sorted([tr_a["region_name"], tr_b["region_name"]]))
            if normalized_pair in comparisons:
                continue

            if tr_a["strand"] == 1:
                exon_5p_a = (tr_a["exon_start"][0], tr_a["exon_stop"][0])
                exon_5p_b = (tr_b["exon_start"][0], tr_b["exon_stop"][0])
            else:
                exon_5p_a = (tr_a["exon_start"][-1], tr_a["exon_stop"][-1])
                exon_5p_b = (tr_b["exon_start"][-1], tr_b["exon_stop"][-1])

            five_prime_overlap = exon_overlap_tup(exon_5p_a, exon_5p_b)
            five_prime_dta, five_prime_ex_3p_dta = gpp.distance_to_agreement(
                exon_5p_a, exon_5p_b
            )
            # distance_2_agreement(
            #     exon_5p_a, exon_5p_b
            # )

            exons_a = [(s, e) for s, e in zip(tr_a["exon_start"], tr_a["exon_stop"])]
            exons_b = [(s, e) for s, e in zip(tr_b["exon_start"], tr_b["exon_stop"])]

            count_90_overlap = gpp.count_overlap(exons_a, exons_b, 0.9)

            # Check if we already calculated similarity for this pair
            type_pair = tuple(sorted([tr_a["so_type"], tr_b["so_type"]]))
            if type_pair in similarity_comparisons:
                type_similarity = similarity_cache[type_pair]
            else:
                type_similarity = rna_type_similarity(
                    tr_a["so_type"], tr_b["so_type"], so_model
                )
                similarity_comparisons.add(type_pair)
                similarity_cache[type_pair] = type_similarity
            comparison = f"{tr_a['region_name']} vs {tr_b['region_name']}"

            if label is not None:
                new_features = pl.DataFrame(
                    {
                        "5p_exon_overlap": five_prime_overlap,
                        "5p_exon_dta": five_prime_dta,
                        "5p_exon_3p_dta": five_prime_ex_3p_dta,
                        "exons_overlapping": count_90_overlap,
                        "strand": int(tr_a["strand"] == tr_b["strand"]),
                        "type_sim": type_similarity,
                        "label": int(label),
                        "comparison": comparison,
                    },
                    schema={
                        "5p_exon_overlap": pl.Float64,
                        "5p_exon_dta": pl.Int64,
                        "5p_exon_3p_dta": pl.Int64,
                        "exons_overlapping": pl.Int64,
                        "strand": pl.Int64,
                        "type_sim": pl.Float64,
                        "label": pl.Int8,
                        "comparison": pl.Utf8,
                    },
                )
            else:
                comparison_features.append(
                    {
                        "5p_exon_overlap": five_prime_overlap,
                        "5p_exon_dta": five_prime_dta,
                        "5p_exon_3p_dta": five_prime_ex_3p_dta,
                        "exons_overlapping": count_90_overlap,
                        "strand": int(tr_a["strand"] == tr_b["strand"]),
                        "type_sim": type_similarity,
                        "comparison": comparison,
                    }
                )
            comparisons.add(normalized_pair)
        new_features = pl.DataFrame(
            comparison_features,
            schema={
                "5p_exon_overlap": pl.Float64,
                "5p_exon_dta": pl.Int64,
                "5p_exon_3p_dta": pl.Int64,
                "exons_overlapping": pl.Int64,
                "strand": pl.Int64,
                "type_sim": pl.Float64,
                "comparison": pl.Utf8,
            },
        )
        features = features.vstack(new_features)
    return features


def polars_preprocessing(transcripts, nearby_distance=1000, label=None):
    chromosomes = sorted(transcripts["chromosome"].unique().to_list())
    all_features = []

    for chrom in (pbar := tqdm(chromosomes[::-1])):
        chrom_transcripts = transcripts.filter(pl.col("chromosome") == chrom)
        pbar.set_description(
            f"Processing chromosome {chrom} ({len(chrom_transcripts)} transcripts)"
        )

        # Skip if too few transcripts to compare
        if len(chrom_transcripts) < 2:
            continue

        features = polars_work_function(chrom_transcripts, nearby_distance)
        if features is not None and len(features) > 0:
            all_features.append(features)

    if all_features:
        all_features_df = pl.concat(all_features).rechunk()
        if label is not None:
            all_features_df = all_features_df.with_columns(pl.lit(label).alias("label"))
        return all_features_df

    return empty_features.clone()


def polars_work_function(transcripts, nearby_distance=1000):
    df = transcripts.with_row_index("idx").sort("region_start")
    ## Minimize memory use by ignoring unneeded columns
    join_cols = [
        "idx",
        "region_start",
        "chromosome",
        "assembly_id",
        "strand",
        "region_name",
        "so_type",
        "exon_start",
        "exon_stop",
    ]
    df_slim = df.select(join_cols).lazy()

    ## Self join with filtering to establish candidate pairs
    ## TODO: need to tweak to get reverse comparison - probably isn't symmetric
    pairs = (
        df_slim.join_where(
            df_slim,
            pl.col("chromosome") == pl.col("chromosome_right"),
            pl.col("assembly_id") == pl.col("assembly_id_right"),
            pl.col("strand") == pl.col("strand_right"),
            pl.col("idx") < pl.col("idx_right"),  # Dedupe: only keep (a,b) not (b,a)
            (pl.col("region_start_right") - pl.col("region_start")).abs()
            <= nearby_distance,
            suffix="_right",
        )
        .with_columns(
            (pl.col("strand") == pl.col("strand_right"))
            .cast(pl.Int8)
            .alias("strand_sim")
        )
        .collect(streaming=True)
    )

    ## Calculate type similarity
    pairs = pairs.with_columns(
        pl.struct(["so_type", "so_type_right"])
        .map_batches(
            lambda s: pl.Series(
                [
                    get_type_similarity(a, b)
                    for a, b in zip(
                        s.struct.field("so_type"), s.struct.field("so_type_right")
                    )
                ]
            ),
            return_dtype=pl.Float32,
        )
        .alias("type_sim")
    )

    ## This replaces the five_prime_overlap = exon_overlap_tup(exon_5p_a, exon_5p_b) call
    pairs = (
        pairs.with_columns(
            [
                # Extract 5' exon coordinates based on strand
                pl.when(pl.col("strand") == 1)
                .then(pl.col("exon_start").list.first())
                .otherwise(pl.col("exon_start").list.last())
                .alias("_5p_start_a"),
                pl.when(pl.col("strand") == 1)
                .then(pl.col("exon_stop").list.first())
                .otherwise(pl.col("exon_stop").list.last())
                .alias("_5p_stop_a"),
                pl.when(pl.col("strand_right") == 1)
                .then(pl.col("exon_start_right").list.first())
                .otherwise(pl.col("exon_start_right").list.last())
                .alias("_5p_start_b"),
                pl.when(pl.col("strand_right") == 1)
                .then(pl.col("exon_stop_right").list.first())
                .otherwise(pl.col("exon_stop_right").list.last())
                .alias("_5p_stop_b"),
            ]
        )
        .with_columns(
            [
                # Normalize coordinates
                pl.min_horizontal("_5p_start_a", "_5p_stop_a").alias("_a_start"),
                pl.max_horizontal("_5p_start_a", "_5p_stop_a").alias("_a_end"),
                pl.min_horizontal("_5p_start_b", "_5p_stop_b").alias("_b_start"),
                pl.max_horizontal("_5p_start_b", "_5p_stop_b").alias("_b_end"),
            ]
        )
        .with_columns(
            [
                (pl.col("_a_end") - pl.col("_a_start")).alias("_len_a"),
                (pl.col("_b_end") - pl.col("_b_start")).alias("_len_b"),
                pl.max_horizontal("_a_start", "_b_start").alias("_overlap_start"),
                pl.min_horizontal("_a_end", "_b_end").alias("_overlap_end"),
            ]
        )
        .with_columns(
            [
                (pl.col("_overlap_end") - pl.col("_overlap_start"))
                .clip(lower_bound=0)
                .alias("_overlap_len"),
            ]
        )
        .with_columns(
            [
                pl.when((pl.col("_len_a") == 0) | (pl.col("_len_b") == 0))
                .then(0.0)
                .otherwise(
                    pl.min_horizontal(
                        pl.col("_overlap_len") / pl.col("_len_a"),
                        pl.col("_overlap_len") / pl.col("_len_b"),
                    )
                )
                .alias("5p_exon_overlap"),
            ]
        )
    ).select(pl.exclude("_*"))

    ## Prepare to calculate other metrics - generate some needed prefeatures
    pairs = pairs.with_columns(
        [
            # Extract 5' exon coordinates based on strand
            pl.when(pl.col("strand") == 1)
            .then(pl.col("exon_start").list.first())
            .otherwise(pl.col("exon_start").list.last())
            .alias("_5p_start_a"),
            pl.when(pl.col("strand") == 1)
            .then(pl.col("exon_stop").list.first())
            .otherwise(pl.col("exon_stop").list.last())
            .alias("_5p_stop_a"),
            pl.when(pl.col("strand_right") == 1)
            .then(pl.col("exon_start_right").list.first())
            .otherwise(pl.col("exon_start_right").list.last())
            .alias("_5p_start_b"),
            pl.when(pl.col("strand_right") == 1)
            .then(pl.col("exon_stop_right").list.first())
            .otherwise(pl.col("exon_stop_right").list.last())
            .alias("_5p_stop_b"),
        ]
    )

    exon_5p_a = list(zip(pairs["_5p_start_a"].to_list(), pairs["_5p_stop_a"].to_list()))
    exon_5p_b = list(zip(pairs["_5p_start_b"].to_list(), pairs["_5p_stop_b"].to_list()))

    # Step 3: Compute distance to agreement (Rust, batched)
    dta, dta_3p = gpp.distance_to_agreement_batch(exon_5p_a, exon_5p_b)
    pairs = pairs.with_columns(
        [
            pl.Series("5p_exon_dta", dta),
            pl.Series("5p_exon_3p_dta", dta_3p),
        ]
    )

    exons_a = [
        list(zip(starts, stops))
        for starts, stops in zip(
            pairs["exon_start"].to_list(), pairs["exon_stop"].to_list()
        )
    ]
    exons_b = [
        list(zip(starts, stops))
        for starts, stops in zip(
            pairs["exon_start_right"].to_list(), pairs["exon_stop_right"].to_list()
        )
    ]
    exon_counts = gpp.count_overlap_batch(exons_a, exons_b, 0.9)
    pairs = pairs.with_columns(pl.Series("exons_overlapping", exon_counts))

    features = pairs.select(
        [
            "5p_exon_overlap",
            "5p_exon_dta",
            "5p_exon_3p_dta",
            "exons_overlapping",
            pl.col("strand_sim").alias("strand"),
            "type_sim",
            pl.concat_str(
                [
                    pl.col("region_name"),
                    pl.lit(" vs "),
                    pl.col("region_name_right"),
                ]
            ).alias("comparison"),
        ]
    )

    return features
