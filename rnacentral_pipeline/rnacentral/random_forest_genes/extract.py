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

from pathlib import Path

import polars as pl
import requests

pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_rows(100)


from tqdm import tqdm


def w_pbar(pbar, func):
    def foo(*args, **kwargs):
        pbar.update(1)
        return func(*args, **kwargs)

    return foo


url = "https://rest.ensembl.org/lookup/id/{0}?content-type=application/json"


def fetch_data(gene):
    r = requests.get(url.format(gene.split(".")[0]))
    data = r.json()

    return {
        "gene_start": data["start"],
        "gene_stop": data["end"],
        "gene_strand": data["strand"],
        "e_chromosome": data["seq_region_name"],
        "assembly": data["assembly_name"],
    }


def locate_genes(df, nearby_distance):
    """
    Loop on each gene, find any genes contained within another gene
    The search genes within 1kb from start or end
    These become the 'nearby' genes to provide the negative examples
    """
    candidates_df = []
    for gene in df.iter_rows(named=True):
        candidates = df.filter(
            pl.col("gene_start").is_between(
                gene["gene_start"] - nearby_distance,
                gene["gene_stop"] + nearby_distance,
            )
            | pl.col("gene_stop").is_between(
                gene["gene_start"] - nearby_distance,
                gene["gene_stop"] + nearby_distance,
            )
            & (pl.col("assembly") == gene["assembly"])
        )
        candidates_df.append(
            {
                "gene": gene["gene"],
                "candidates": [gene["gene"]],
                "labels": [1],
                "contained": [True],
            }
        )
        for c in candidates.iter_rows(named=True):
            if c["gene"] == gene["gene"]:
                continue
            candidates_df[-1]["candidates"].append(c["gene"])
            candidates_df[-1]["labels"].append(0)
            contained = (
                c["gene_start"] > gene["gene_start"]
                and c["gene_stop"] < gene["gene_stop"]
            )
            candidates_df[-1]["contained"].append(contained)

    candidates_df = pl.DataFrame(candidates_df)
    return candidates_df


def genes_transcripts(
    dataframe,
    taxid,
    output_transcripts,
    output_genes,
    ensembl_only,
    pre_gene_set,
    nearby_distance,
):
    data = pl.scan_parquet(dataframe).filter(pl.col("taxid") == taxid)

    if ensembl_only:
        data = data.filter(pl.col("gene").str.starts_with("ENS"))
    data = data.group_by(
        ["assembly_id", "gene", "region_name", "strand"], maintain_order=True
    ).agg(
        [
            pl.col("urs_taxid").first(),
            pl.col("chromosome").first(),
            # pl.col("strand"),
            pl.col("region_start").first(),
            pl.col("region_stop").first(),
            pl.col("exon_start"),
            pl.col("exon_stop"),
            pl.col("so_rna_type").first().alias("so_type"),
        ]
    )

    if pre_gene_set is not None:
        if Path(pre_gene_set).exists():
            genes = pl.read_parquet(pre_gene_set)
        else:
            genes = data.select(pl.col("gene").unique()).collect()
            if genes.height == 0:
                print("No genes found with selected parameters, try something else")
                exit()
            print(genes)
            pbar = tqdm(total=len(genes), desc="fetching gene coords", colour="green")
            genes = genes.with_columns(
                res=pl.col("gene").map_elements(w_pbar(pbar, fetch_data))
            ).unnest("res")
            genes.write_parquet(pre_gene_set)
            pbar.close()
    else:
        print("You need to download the gene set. Give a filepath to --pre_gene_set")

    print("Calculating nearby genes")
    genes = genes.filter(pl.col("gene_start").is_not_null())
    candidates = genes.group_by(["assembly", "e_chromosome"]).map_groups(
        lambda x: locate_genes(x, nearby_distance)
    )

    if output_genes is not None:
        ## Because we ask to maintain order in the groupby, we can't stream
        ## That's unfortunate, but maintaining order is needed to be sure the
        ## exons are in the right order
        candidates.write_parquet(output_genes)

    if output_transcripts is not None:
        data.collect().write_parquet(output_transcripts)
