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

import bisect
import multiprocessing as mp
import time
from pathlib import Path

## Use rust preprocessing code
## Must be built using the Makefile in an activated environment
import gene_preprocessing as gpp
import numpy as np
import polars as pl
from gensim.models import Word2Vec
from line_profiler import profile
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
        "type_sim": pl.Float64,
        "comparison": pl.Utf8,
    },
)


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


@profile
def compare_transcripts(transcripts_a, transcripts_b, so_model, label=0):
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


@profile
def identify_nearby_transcripts(transcripts, so_model, nearby_distance=1000):
    """
    Find transcripts within 1kb of each other to calculate feature sets for and calculate features

    Args:
        transcripts (pl.DataFrame): DataFrame of transcripts
        nearby_distance (int): Distance to search for nearby transcripts
    Returns:
        pl.DataFrame: DataFrame of features
    """

    def process_transcripts():
        for transcript_a in tqdm(
            transcripts.iter_rows(named=True),
            total=transcripts.height,
            desc="Calculating features...",
        ):
            candidates = transcripts.filter(
                (
                    pl.col("region_start").is_between(
                        transcript_a["region_start"] - nearby_distance,
                        transcript_a["region_stop"] + nearby_distance,
                    )
                    | pl.col("region_stop").is_between(
                        transcript_a["region_start"] - nearby_distance,
                        transcript_a["region_stop"] + nearby_distance,
                    )
                )
                & (pl.col("region_name") != transcript_a["region_name"])
                & (pl.col("chromosome") == transcript_a["chromosome"])
                & (pl.col("assembly_id") == transcript_a["assembly_id"])
            )

            if len(candidates) > 0:
                yield from compare_transcripts(
                    pl.DataFrame([transcript_a]), candidates, so_model, label=None
                )

    return pl.DataFrame(process_transcripts())


def identify_nearby_transcripts_sorted(transcripts, so_model, nearby_distance=1000):
    """Optimized version using sorting and binary search"""

    sorted_keys = (
        transcripts.group_by(["chromosome", "assembly_id", "strand"])
        .len()
        .sort("len", descending=True)
        .select(["chromosome", "assembly_id", "strand"])
    )
    # print(sorted_keys)
    # # Group by chromosome and assembly_id
    # grouped = transcripts.group_by(["chromosome", "assembly_id", "strand"])
    # print(grouped.df)
    # exit()
    ## This should not have labels in because we're in the inference side
    features = pl.DataFrame(
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
            "type_sim": pl.Float64,
            "comparison": pl.Utf8,
        },
    )

    for row in sorted_keys.iter_rows(named=True):
        group_df = transcripts.filter(
            (pl.col("chromosome") == row["chromosome"])
            & (pl.col("assembly_id") == row["assembly_id"])
            & (pl.col("strand") == row["strand"])
        )
        # Sort by region_start for binary search
        sorted_df = group_df.sort("region_start")
        sorted_data = sorted_df.to_dicts()

        group_key = f"{row['assembly_id']}::{row['chromosome']}::{row['strand']}"

        for i, transcript_a in enumerate(
            tqdm(sorted_data, desc=f"Processing {group_key}")
        ):

            # Use binary search to find the range of potentially nearby transcripts
            search_start = transcript_a["region_start"] - nearby_distance
            search_end = transcript_a["region_start"] + nearby_distance

            # Find candidates within the search range
            candidates_data = []
            for j, transcript_b in enumerate(sorted_data):
                if i == j:  # Skip self
                    continue

                # Early termination if we've passed the search range
                if transcript_b["region_start"] > search_end:
                    break

                # Check if within range
                if (
                    transcript_b["region_start"] <= search_end
                    and transcript_b["region_start"] >= search_start
                ):
                    candidates_data.append(transcript_b)

            if candidates_data:
                candidates = pl.DataFrame(candidates_data)
                f = compare_transcripts(
                    pl.DataFrame([transcript_a]), candidates, so_model, label=None
                )
                features = features.vstack(f)
    return features.rechunk()


def _process_group_worker(args):
    """
    Worker function for processing a single chromosome/assembly group.

    Args:
        args: Tuple containing (group_key, group_data, so_model_path, nearby_distance)

    Returns:
        List of feature dictionaries for the group
    """
    group_key, group_data, so_model_path, nearby_distance = args

    # Load the model in each worker process
    so_model = Word2Vec.load(so_model_path)

    # Convert group data back to DataFrame
    group_df = pl.DataFrame(group_data)

    # Sort by region_start for binary search
    sorted_df = group_df.sort("region_start")
    sorted_data = sorted_df.to_dicts()

    features = []
    print(f"starting group {group_key} with {len(sorted_data)} transcripts")

    for i, transcript_a in enumerate(sorted_data):
        # Use binary search to find the range of potentially nearby transcripts
        search_start = transcript_a["region_start"] - nearby_distance
        search_end = transcript_a["region_stop"] + nearby_distance

        # Find candidates within the search range
        candidates_data = []
        for j, transcript_b in enumerate(sorted_data):
            if i == j:  # Skip self
                continue

            # Early termination if we've passed the search range
            if transcript_b["region_start"] > search_end:
                break

            # Check if within range
            if (
                transcript_b["region_start"] <= search_end
                and transcript_b["region_stop"] >= search_start
            ):
                candidates_data.append(transcript_b)

        if candidates_data:
            candidates = pl.DataFrame(candidates_data)
            f = compare_transcripts(
                pl.DataFrame([transcript_a]), candidates, so_model, label=None
            )
            features.extend(f)
    print(f"Processed group {group_key} in {time.time() - start:.2f} seconds")
    return features


def _process_group_worker_optimized(args):
    """
    Optimized worker function for processing a single chromosome/assembly group.
    Fixes O(n²) search and expensive Cartesian product issues.

    Args:
        args: Tuple containing (group_key, group_data, so_model_path, nearby_distance)

    Returns:
        List of feature dictionaries for the group
    """
    group_key, group_data, so_model, nearby_distance = args

    # Load the model in each worker process
    # so_model = Word2Vec.load(so_model_path)

    # Convert group data back to DataFrame
    group_df = pl.DataFrame(group_data)

    # Sort by region_start for binary search
    sorted_df = group_df.sort("region_start")
    sorted_data = sorted_df.to_dicts()

    # Pre-extract region starts for binary search optimization
    region_starts = sorted_df.get_column("region_start").to_list()

    features = []
    print(f"starting group {group_key} with {len(sorted_data)} transcripts")
    if group_key[0] == "MT":
        print(
            sorted_df.with_columns(exon_count=pl.col("exon_stop").list.len()).describe()
        )
    start = time.time()
    for i, transcript_a in enumerate(sorted_data):
        # Use proper binary search to find the range of potentially nearby transcripts
        search_start = transcript_a["region_start"] - nearby_distance
        search_end = transcript_a["region_stop"] + nearby_distance

        # Find the range using binary search - much faster than O(n) iteration
        left_idx = bisect.bisect_left(region_starts, search_start)
        right_idx = bisect.bisect_right(region_starts, search_end)

        # Find candidates within the search range
        candidates_data = []
        for j in range(left_idx, min(right_idx, len(sorted_data))):
            if i == j:  # Skip self
                continue

            transcript_b = sorted_data[j]
            # Check if within range (region_stop check still needed)
            if transcript_b["region_stop"] >= search_start:
                candidates_data.append(transcript_b)

        if candidates_data:
            candidates = pl.DataFrame(candidates_data)
            f = compare_transcripts(
                pl.DataFrame([transcript_a]), candidates, so_model, label=None
            )
            features.extend(f)
    print(f"Processed group {group_key} in {time.time() - start:.2f} seconds")

    return features


def identify_nearby_transcripts_sorted_parallel(
    transcripts, so_model_path, nearby_distance=1000, n_processes=None
):
    """
    Parallelized version using multiprocessing across chromosome/assembly groups.

    Args:
        transcripts (pl.DataFrame): DataFrame of transcripts
        so_model_path (str): Path to the Word2Vec model file
        nearby_distance (int): Distance to search for nearby transcripts
        n_processes (int, optional): Number of processes to use. Defaults to CPU count.

    Returns:
        pl.DataFrame: DataFrame of features
    """
    if n_processes is None:
        n_processes = min(mp.cpu_count(), 8)  # Cap at 8 to avoid excessive memory usage

    # Group by chromosome and assembly_id
    grouped = transcripts.group_by(["chromosome", "assembly_id"])

    # Prepare arguments for worker processes
    # Convert DataFrames to dictionaries to avoid pickling issues
    worker_args = []
    so_model = Word2Vec.load(so_model_path)
    for group_key, group_df in grouped:
        group_data = group_df.to_dicts()
        if len(group_data) < 1000:
            # Skip groups with less than 2 transcripts
            continue
        worker_args.append((group_key, group_data, so_model, nearby_distance))

    # Sort by group size (largest first) for better load balancing
    worker_args.sort(key=lambda x: len(x[1]), reverse=True)

    all_features = []

    if len(worker_args) == 1 or n_processes == 1:
        # Fall back to sequential processing for single group or single process
        print("Processing sequentially (single group or single process)")
        for args in tqdm(worker_args, desc="Processing groups"):
            features = _process_group_worker_optimized(args)
            all_features.extend(features)
    else:
        # Use multiprocessing for multiple groups
        print(f"Processing {len(worker_args)} groups using {n_processes} processes")

        with mp.Pool(processes=n_processes) as pool:
            # Use imap for progress tracking
            results = list(
                tqdm(
                    pool.map(_process_group_worker_optimized, worker_args, chunksize=1),
                    total=len(worker_args),
                    desc="Processing groups",
                )
            )
            # Flatten results
            for features in results:
                all_features.extend(features)

    return pl.DataFrame(all_features)


def run_preprocessing(
    transcripts_file,
    regions_data,
    so_model_path,
    nearby_distance,
    use_parallel=True,
    n_processes=None,
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
            transcripts = data.add_assembly_region_ids(transcripts, regions_data)

        if use_parallel:
            features = identify_nearby_transcripts_sorted_parallel(
                transcripts, so_model_path, nearby_distance, n_processes
            )
        else:
            # Load model for sequential processing
            so_model = Word2Vec.load(so_model_path)
            so_vec_normalised = {
                key: so_model.wv[key] / np.linalg.norm(so_model.wv[key])
                for key in so_model.wv.key_to_index
            }
            features = identify_nearby_transcripts_sorted(
                transcripts, so_vec_normalised, nearby_distance
            )
    else:
        features = pl.DataFrame()

    return features


if __name__ == "__main__":
    transcripts = pl.read_parquet("/Users/agreen/code/rnc-genes/human_transcripts.pq")

    so_model = Word2Vec.load(
        "/Users/agreen/code/rnc-genes/node2vec/so_embedding_model.emb"
    )

    identify_nearby_transcripts_sorted(transcripts, so_model)
