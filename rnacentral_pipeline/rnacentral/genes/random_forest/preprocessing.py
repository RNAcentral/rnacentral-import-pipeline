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


import tempfile
from functools import lru_cache
from pathlib import Path

## Use rust preprocessing code
## Must be built using the Makefile in an activated environment
import gene_preprocessing as gpp
import numpy as np
import polars as pl
from gensim.models import Word2Vec
from tqdm import tqdm
import logging
LOGGER = logging.getLogger(__name__)

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
        "strand": pl.Int8,
        "type_sim": pl.Float32,
        "comparison": pl.Utf8,
    },
)

obsolete_so_terms = {
    "SO:0001171": "SO:0002345"
}

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
    if type_a in obsolete_so_terms:
        LOGGER.warning(f"Obsolete SO term {type_a} replaced with {obsolete_so_terms[type_a]}")
        type_a = obsolete_so_terms[type_a]
    if type_b in obsolete_so_terms:
        LOGGER.warning(f"Obsolete SO term {type_b} replaced with {obsolete_so_terms[type_b]}")
        type_b = obsolete_so_terms[type_b]
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
    output_path=None,
):
    """
    Generate pairwise comparison features from a transcripts parquet file.

    When ``output_path`` is given the features are streamed straight to that
    parquet file (peak memory bounded by a single chromosome) and the number of
    rows written is returned. Otherwise the full feature DataFrame is returned.
    """
    from rnacentral_pipeline.rnacentral.genes.random_forest import data

    transcripts = pl.read_parquet(transcripts_file)
    ## Filter out piRNAs
    transcripts = transcripts.filter(pl.col("so_type") != "SO:0001035")

    if transcripts.height == 0:
        features = empty_features.clone()
        if output_path is not None:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            features.write_parquet(output_path)
            return 0
        return features

    # Check if region_ids are present, fetch if needed
    if "region_id" not in transcripts.columns:
        if not Path(regions_data).exists():
            raise ValueError(
                "Region IDs not found in transcripts file and no database connection provided. "
                "Please provide --conn_str to fetch region IDs from database."
            )
        transcripts = data.add_assembly_region_ids(transcripts, regions_data).select(
            pl.exclude("assembly_id_right")
        )

    init_so_model(so_model_path)

    if output_path is not None:
        return sink_preprocessing(transcripts, output_path, nearby_distance)

    return polars_preprocessing(transcripts, nearby_distance)


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
                "strand": pl.Int8,
                "type_sim": pl.Float32,
                "comparison": pl.Utf8,
            },
        )
        features = features.vstack(new_features)
    return features


def _spill_chromosome_features(transcripts, tmpdir, nearby_distance):
    """
    Run polars_work_function per chromosome, writing each chromosome's feature
    frame to its own parquet file under ``tmpdir``.

    Returns the list of written file paths. Each chromosome's features are
    *streamed* straight to its parquet file via the streaming engine, so the
    full pair table for a chromosome is never held in memory -- peak RAM is
    bounded by the streaming engine's batch size, not the chromosome size.
    """
    chromosomes = sorted(transcripts["chromosome"].unique().to_list())
    paths = []

    for chrom in (pbar := tqdm(chromosomes[::-1])):
        chrom_transcripts = transcripts.filter(pl.col("chromosome") == chrom)
        pbar.set_description(
            f"Processing chromosome {chrom} ({chrom_transcripts.height} transcripts)"
        )

        # Skip if too few transcripts to compare
        if chrom_transcripts.height < 2:
            continue

        path = tmpdir / f"chrom_{len(paths):05d}.parquet"
        _work_plan(chrom_transcripts, nearby_distance).sink_parquet(
            path, engine="streaming"
        )
        paths.append(path)
        ## Drop the per-chromosome transcripts before the next (potentially
        ## large) chromosome is processed.
        del chrom_transcripts

    return paths


def polars_preprocessing(transcripts, nearby_distance=1000, label=None):
    """
    Build the full feature table in memory.

    Each chromosome is streamed to a temporary parquet file (bounded memory),
    then the per-chromosome files are read back and concatenated. Only the final
    concatenated frame is held in full.

    For the production path prefer ``sink_preprocessing``, which streams the
    result straight to disk and never materialises the whole table.
    """
    with tempfile.TemporaryDirectory(prefix="gene_preproc_", dir=".") as tmp:
        paths = _spill_chromosome_features(transcripts, Path(tmp), nearby_distance)

        if not paths:
            features = empty_features.clone()
        else:
            features = pl.read_parquet([str(p) for p in paths]).rechunk()

        if label is not None:
            features = features.with_columns(pl.lit(label).alias("label"))

        return features


def sink_preprocessing(transcripts, output_path, nearby_distance=1000, label=None):
    """
    Like ``polars_preprocessing`` but stream the result straight to
    ``output_path`` (a parquet file) instead of returning it. Returns the number
    of feature rows written.

    Each chromosome is itself streamed to a part file (see
    ``_spill_chromosome_features``); the parts are then combined with a lazy
    scan + sink, so the full table is never materialised in memory.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ## Spill next to the output file rather than in the system temp dir: on the
    ## cluster /tmp is a small node-local filesystem that fills up, and the
    ## output directory is on the same (large) filesystem as the final result.
    with tempfile.TemporaryDirectory(
        prefix="gene_preproc_", dir=output_path.parent
    ) as tmp:
        paths = _spill_chromosome_features(transcripts, Path(tmp), nearby_distance)

        if not paths:
            out = empty_features.clone()
            if label is not None:
                out = out.with_columns(pl.lit(label).alias("label"))
            out.write_parquet(output_path)
            return 0

        lazy = pl.scan_parquet([str(p) for p in paths])
        if label is not None:
            lazy = lazy.with_columns(pl.lit(label).alias("label"))
        lazy.sink_parquet(output_path)

    ## Row count without reading the data back into memory.
    return pl.scan_parquet(output_path).select(pl.len()).collect().item()


def _type_similarity_lookup(transcripts):
    """
    Precompute a (so_type, so_type_right) -> type_sim lookup as a LazyFrame.

    The set of SO types in a chromosome is tiny (tens), so the full cross
    product is small. Joining this into the plan replaces the per-row Python UDF
    that used to compute type similarity, which is what forced the whole pair
    table to be materialised. ``get_type_similarity`` is symmetric, so the value
    is the same regardless of which side a type appears on.
    """
    types = transcripts["so_type"].unique().to_list()
    return pl.DataFrame(
        [(a, b, float(get_type_similarity(a, b))) for a in types for b in types],
        schema={
            "so_type": pl.Utf8,
            "so_type_right": pl.Utf8,
            "type_sim": pl.Float32,
        },
        orient="row",
    ).lazy()


def _work_plan(transcripts, nearby_distance=1000):
    """
    Build the per-chromosome feature computation as a single lazy plan.

    Every feature is expressed natively in polars (no Python UDFs, no Rust
    round-trips), so the plan can be run under the streaming engine and sunk
    straight to disk -- peak memory is bounded by the streaming engine's batch
    size rather than the full pair table. ``polars_work_function`` collects this
    plan; ``_spill_chromosome_features`` sinks it.

    The exon-overlap count is computed natively here and has been fuzz-tested
    against the Rust ``count_overlap_batch`` (including exact-0.9 boundaries) for
    exact equivalence.
    """
    ## A stable index in region_start order defines the (a, b) orientation used
    ## below: idx < idx_right keeps only one of each pair, with the lower
    ## region_start as 'a'.
    df = transcripts.with_row_index("idx").sort("region_start").with_columns(
        ## 5' exon coordinates: first exon on the + strand, last on the -.
        fivep_start=pl.when(pl.col("strand") == 1)
        .then(pl.col("exon_start").list.first())
        .otherwise(pl.col("exon_start").list.last()),
        fivep_stop=pl.when(pl.col("strand") == 1)
        .then(pl.col("exon_stop").list.first())
        .otherwise(pl.col("exon_stop").list.last()),
    )

    type_lut = _type_similarity_lookup(transcripts)

    join_cols = [
        "idx",
        "region_start",
        "chromosome",
        "assembly_id",
        "strand",
        "region_name",
        "so_type",
        "fivep_start",
        "fivep_stop",
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
        ## strand_sim is always 1 here (the join requires strand == strand_right)
        ## but is computed the same way the row-wise reference does, for parity.
        .with_columns(
            (pl.col("strand") == pl.col("strand_right"))
            .cast(pl.Int8)
            .alias("strand_sim")
        )
        ## 5' exon overlap and distance-to-agreement from the scalar 5' coords.
        ## distance_to_agreement is just the absolute difference of the
        ## (start, stop) pair.
        .with_columns(
            [
                pl.min_horizontal("fivep_start", "fivep_stop").alias("_a_start"),
                pl.max_horizontal("fivep_start", "fivep_stop").alias("_a_end"),
                pl.min_horizontal("fivep_start_right", "fivep_stop_right").alias(
                    "_b_start"
                ),
                pl.max_horizontal("fivep_start_right", "fivep_stop_right").alias(
                    "_b_end"
                ),
                (pl.col("fivep_start") - pl.col("fivep_start_right"))
                .abs()
                .alias("5p_exon_dta"),
                (pl.col("fivep_stop") - pl.col("fivep_stop_right"))
                .abs()
                .alias("5p_exon_3p_dta"),
            ]
        )
        .with_columns(
            [
                (pl.col("_a_end") - pl.col("_a_start")).alias("_len_a"),
                (pl.col("_b_end") - pl.col("_b_start")).alias("_len_b"),
                (
                    pl.min_horizontal("_a_end", "_b_end")
                    - pl.max_horizontal("_a_start", "_b_start")
                )
                .clip(lower_bound=0)
                .alias("_overlap_len"),
            ]
        )
        .with_columns(
            pl.when((pl.col("_len_a") == 0) | (pl.col("_len_b") == 0))
            .then(0.0)
            .otherwise(
                pl.min_horizontal(
                    pl.col("_overlap_len") / pl.col("_len_a"),
                    pl.col("_overlap_len") / pl.col("_len_b"),
                )
            )
            .alias("5p_exon_overlap")
        )
    )
    pairs = (
        pairs.with_columns(
            _n=pl.min_horizontal(
                pl.col("exon_start").list.len(),
                pl.col("exon_start_right").list.len(),
            )
        )
        .with_columns(
            _eas=pl.col("exon_start").list.slice(0, pl.col("_n")),
            _eae=pl.col("exon_stop").list.slice(0, pl.col("_n")),
            _ebs=pl.col("exon_start_right").list.slice(0, pl.col("_n")),
            _ebe=pl.col("exon_stop_right").list.slice(0, pl.col("_n")),
        )
        .with_columns(
            ## per-exon overlap start = max(_eas, _ebs), end = min(_eae, _ebe)
            _ov_start=pl.col("_eas")
            + (pl.col("_ebs") - pl.col("_eas")).list.eval(
                pl.element().clip(lower_bound=0)
            ),
            _ov_end=pl.col("_eae")
            - (pl.col("_eae") - pl.col("_ebe")).list.eval(
                pl.element().clip(lower_bound=0)
            ),
            _len_ea=pl.col("_eae") - pl.col("_eas"),
            _len_eb=pl.col("_ebe") - pl.col("_ebs"),
        )
        .with_columns(
            _ov_len=(pl.col("_ov_end") - pl.col("_ov_start")).list.eval(
                pl.element().clip(lower_bound=0)
            ),
        )
        .with_columns(
            ## ratios computed in Float32 to match the Rust implementation.
            _ratio_a=pl.col("_ov_len").cast(pl.List(pl.Float32))
            / pl.col("_len_ea").cast(pl.List(pl.Float32)),
            _ratio_b=pl.col("_ov_len").cast(pl.List(pl.Float32))
            / pl.col("_len_eb").cast(pl.List(pl.Float32)),
        )
        .with_columns(
            _min_ratio=pl.col("_ratio_a")
            - (pl.col("_ratio_a") - pl.col("_ratio_b")).list.eval(
                pl.element().clip(lower_bound=0)
            )
        )
        .with_columns(
            ## count positions with a real overlap (>0) whose smaller ratio >= 0.9
            exons_overlapping=(
                pl.col("_ov_len").list.eval(pl.element() > 0).cast(pl.List(pl.Int8))
                * pl.col("_min_ratio")
                .list.eval(pl.element() >= pl.lit(0.9, dtype=pl.Float32))
                .cast(pl.List(pl.Int8))
            )
            .list.sum()
            .cast(pl.Int64)
        )
    )

    ## type similarity via the precomputed lookup join (replaces a Python UDF).
    pairs = pairs.join(type_lut, on=["so_type", "so_type_right"], how="left")

    return pairs.select(
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


def polars_work_function(transcripts, nearby_distance=1000):
    """
    Compute the feature table for a single chromosome's transcripts.

    Thin wrapper that collects the lazy plan from ``_work_plan`` under the
    streaming engine. ``_spill_chromosome_features`` sinks the same plan instead
    of collecting it.
    """
    return _work_plan(transcripts, nearby_distance).collect(engine="streaming")
