# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

from itertools import product

import numpy as np
import polars as pl
from gensim.models import Word2Vec

pl.Config.set_tbl_cols(-1)


so_model = Word2Vec.load("/Users/agreen/code/node2vec/so_embedding_model.emb")


def nearby(vals, distance=1e3):
    """
    Require sorted dataframe
    look at this entry and the 3 before it, are they within a given distance to become 'nearby'
    Should consider start + stop, strand and chromosome, nothing else
    try 1kb as the distance?
    """
    if vals[0] < distance:
        return True
    else:
        return False


def exon_overlap_tup(exon_a, exon_b):
    return exon_overlap(exon_a[0], exon_a[1], exon_b[0], exon_b[1])


def distance_2_agreement(exon_a, exon_b):
    """
    Only look at the squared distance between the start coordinates
    """
    dist_5p = np.sqrt((exon_a[0] - exon_b[0]) ** 2)
    dist_3p = np.sqrt((exon_a[1] - exon_b[1]) ** 2)
    return (dist_5p, dist_3p)


def count_overlapping(starts_a, stops_a, starts_b, stops_b):
    """
    Count the number of exons with
    """


def exon_overlap(exon_a_start, exon_a_end, exon_b_start, exon_b_end):
    """
    Calculates overlap of exon b with exon a
    """
    overlap_start = max(exon_a_start, exon_b_start)
    overlap_stop = min(exon_a_end, exon_b_end)
    overlap_length = max(0, overlap_stop - overlap_start)

    length_a = exon_a_end - exon_a_start
    length_b = exon_b_end - exon_b_start

    if length_a > 0:
        overlap_a = overlap_length / length_a
    else:
        overlap_a = 999

    if length_b > 0:
        overlap_b = overlap_length / length_b
    else:
        overlap_b = 999
    overlap = min(overlap_a, overlap_b)

    if overlap > 1:
        overlap = 0

    return overlap


def rna_type_similarity(so_type_a, so_type_b):
    """
    Compute dot product similarity between transcript type vectors
    """
    so_vec_a = so_model.wv[so_type_a]
    so_vec_b = so_model.wv[so_type_b]

    sim = np.dot(so_vec_a, so_vec_b) / (
        np.linalg.norm(so_vec_a) * np.linalg.norm(so_vec_b)
    )
    return sim


def compare_transcripts(transcripts_a, transcripts_b, label=0):
    comparisons = []
    features = []
    for tr_a in transcripts_a.iter_rows(named=True):
        for tr_b in transcripts_b.iter_rows(named=True):
            if tr_a["region_name"] == tr_b["region_name"]:
                continue
            if (tr_a["region_name"], tr_b["region_name"]) in comparisons or (
                tr_b["region_name"],
                tr_a["region_name"],
            ) in comparisons:
                continue

            if tr_a["strand"] == 1:
                exon_5p_a = (tr_a["exon_start"][0], tr_a["exon_stop"][0])
            else:
                exon_5p_a = (tr_a["exon_start"][-1], tr_a["exon_stop"][-1])

            if tr_b["strand"] == 1:
                exon_5p_b = (tr_b["exon_start"][0], tr_b["exon_stop"][0])
            else:
                exon_5p_b = (tr_b["exon_start"][-1], tr_b["exon_stop"][-1])

            five_prime_overlap = exon_overlap_tup(exon_5p_a, exon_5p_b)
            five_prime_dta, five_prime_ex_3p_dta = distance_2_agreement(
                exon_5p_a, exon_5p_b
            )

            exons_a = [(s, e) for s, e in zip(tr_a["exon_start"], tr_a["exon_stop"])]
            exons_b = [(s, e) for s, e in zip(tr_b["exon_start"], tr_b["exon_stop"])]
            count_90_overlap = int(
                sum(
                    [
                        exon_overlap_tup(xa, xb) > 0.9
                        for xa, xb in product(exons_a, exons_b)
                    ]
                )
            )

            type_similarity = rna_type_similarity(tr_a["so_type"], tr_b["so_type"])
            comparison = f"{tr_a['region_name']} vs {tr_b['region_name']}"

            if label is not None:
                features.append(
                    {
                        "5p_exon_overlap": five_prime_overlap,
                        "5p_exon_dta": five_prime_dta,
                        "5p_exon_3p_dta": five_prime_ex_3p_dta,
                        "exons_overlapping": count_90_overlap,
                        "strand": int(tr_a["strand"] == tr_b["strand"]),
                        "type_sim": type_similarity,
                        "label": int(label),
                        "comparison": comparison,
                    }
                )
            else:
                features.append(
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
            comparisons.extend(
                [
                    (tr_a["region_name"], tr_b["region_name"]),
                    (tr_b["region_name"], tr_a["region_name"]),
                ]
            )
    return features


def transcripts(candidates, transcript_data, feature_output):
    candidates = pl.read_parquet(candidates)
    t_data = pl.read_parquet(
        transcript_data
    )  # .with_columns(pl.col("so_type").list.first())

    features = []
    for row in candidates.iter_rows(named=True):
        ## The extract genes step has already guaranteed that genes and candidates are from the same assembly
        cand_genes = row["candidates"]
        cand_labels = row["labels"]
        row_transcripts = t_data.filter(pl.col("gene") == row["gene"])

        for cand, label in zip(cand_genes, cand_labels):
            transcripts = t_data.filter(pl.col("gene") == cand)

            ## Mark up all transcripts in this gene as being in the gene
            features.extend(
                compare_transcripts(transcripts, row_transcripts, label=label)
            )

    example_data = pl.DataFrame(features).unique("comparison")

    example_data.write_parquet(feature_output)
