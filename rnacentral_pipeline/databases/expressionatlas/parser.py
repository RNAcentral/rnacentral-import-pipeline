# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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
import json
import operator as op
import pathlib
import typing as ty

import polars as pl

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.expressionatlas import sdrf

from . import helpers


def parse_differential(analytics, sdrf_path, lookup):
    """
    Join the analytics against the lookup to get only rows for genes we care about

    Analytics is filtered by:
    - non null p-value. Note the non-strict cast to coerce NA into null
    - log2foldchange >= 1
    """
    analytics_data = pl.read_csv(
        analytics,
        separator="\t",
    ).with_columns(pl.selectors.contains("p-value").cast(pl.Float32, strict=False))
    if analytics_data.height == 0:
        raise ValueError(f"Analytics data for {analytics} was empty, abort parsing")
    analytics_data = (
        analytics_data.filter(
            pl.any_horizontal(pl.selectors.contains("p-value").is_not_null())
        )
        .filter(pl.any_horizontal(pl.selectors.contains("p-value").lt(0.05)))
        .filter(pl.any_horizontal(pl.selectors.contains("log2foldchange").abs().ge(1)))
        .with_columns(experiment=pl.lit(analytics.parents[-1]))
    )
    if "Gene ID" in analytics_data.columns:
        analytics_data = analytics_data.rename({"Gene ID": "GeneID"})

    analytics_data = analytics_data.select(["GeneID", "experiment"])

    sdrf_data = sdrf.parse_condensed_sdrf(sdrf_path)

    organisms = (
        sdrf_data.filter(pl.col("feat_type") == "organism").select("ontology").unique()
    )
    taxids = organisms.with_columns(
        taxid=pl.col("ontology").str.split("NCBITaxon_").list.last().cast(pl.Int64)
    ).select("taxid")

    ## The skip_rows=2 here is because psql will write the CREATE and COPY acknowledgements in front of the csv data
    lookup_data = (
        pl.scan_csv(lookup, skip_rows=2)
        .unique(["urs_taxid", "gene"])
        .select(["urs_taxid", "taxid", "gene"])
        .join(taxids.lazy(), on="taxid")
    )
    output_data = (
        lookup_data.join(analytics_data.lazy(), left_on="gene", right_on="GeneID")
        .select(["urs_taxid", "experiment"])
        .collect(streaming=True)
    )

    return output_data


def parse_baseline(tpms, sdrf_path, lookup):
    """
    Parse and filter the baseline tpms file, then get the genes we care about

    - tpms is parsed to get to a median tpm from a lit of replicates (like 1,2,3,4,5 )
    - Then we filter to only keep things with median > 0.5
    - Then join against the lookup table for things with this taxid to extract gene-experiment-urs linkage
    - Return dataframe with only necessary fields

    """

    tpms_data = (
        pl.read_csv(tpms, separator="\t")
        .with_columns(
            pl.selectors.matches("^g\d+")
            .str.split(",")
            .list.eval(pl.element().cast(pl.Float64))
        )
        .with_columns(pl.selectors.matches("^g\d+").list.median())
        .filter(pl.any_horizontal(pl.selectors.matches("^g\d+").gt(0.5)))
        .with_columns(experiment=pl.lit(tpms.parents[-1]))
    )
    if "Gene ID" in tpms_data.columns:
        tpms_data = tpms_data.rename({"Gene ID": "GeneID"})

    tpms_data = tpms_data.select(["GeneID", "experiment"])

    sdrf_data = sdrf.parse_condensed_sdrf(sdrf_path)

    organisms = (
        sdrf_data.filter(pl.col("feat_type") == "organism").select("ontology").unique()
    )
    taxids = organisms.with_columns(
        taxid=pl.col("ontology").str.split("NCBITaxon_").list.last().cast(pl.Int64)
    ).select("taxid")

    lookup_data = (
        pl.scan_csv(lookup, skip_rows=2)
        .unique(["urs_taxid", "gene"])
        .select(["urs_taxid", "taxid", "gene"])
        .join(taxids.lazy(), on="taxid")
    )
    output_data = (
        lookup_data.join(tpms_data.lazy(), left_on="gene", right_on="GeneID")
        .select(["urs_taxid", "experiment"])
        .collect(streaming=True)
    )

    return output_data


def parse(handle, lookup):
    """
    Process the jsonlines output from Rust into entries.

    The jsonlines output is already grouped by geneID and urs taxid so this
    should give us the transcript level linkage we're after without any further
    processing.
    """
    ## Group on urs to get a list of matching experiments
    grouped_data = (
        pl.scan_ndjson(handle).group_by("urs_taxid").agg(pl.col("experiment"))
    )

    lookup_data = pl.scan_csv(lookup, skip_rows=2)
    ## Join against lookup to get the rest of the required information
    hits = grouped_data.join(lookup_data, on="urs_taxid", how="inner").collect(
        streaming=True
    )
    for hit in hits.iter_rows(named=True):
        for experiment in hit["experiment"]:
            yield helpers.as_entry(hit, experiment)
