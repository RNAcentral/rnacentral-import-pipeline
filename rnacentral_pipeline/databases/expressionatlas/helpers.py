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

import logging
from pathlib import Path

import polars as pl

from rnacentral_pipeline.databases.data import Entry, Exon, SequenceRegion
from rnacentral_pipeline.databases.expressionatlas import sdrf
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs

LOGGER = logging.getLogger(__name__)


def accession(info):
    return "EXPRESSIONATLAS:" + info["gene"]


def primary_id(info):
    return "EXPRESSIONATLAS:" + info["gene"]


def taxid(info):
    taxid = info["taxid"]
    return int(taxid)


def species(info):
    return phy.species(info["taxid"])


def lineage(info):
    return phy.lineage(info["taxid"])


def common_name(info):
    return phy.common_name(info["taxid"])


def url(info):
    return "https://www.ebi.ac.uk/gxa/genes/" + info["gene"]


def region_builder(info):
    return []


def references(interactions):
    refs = set()
    for interaction in interactions:
        refs.update(interaction.publications)
    refs.add(pubs.reference(24234451))
    return list(refs)


def rna_type(type_str):
    if type_str is None:
        return "SO:0000655"
    else:
        return type_str


def as_entry(info):
    return Entry(
        primary_id=primary_id(info),
        accession=accession(info),
        ncbi_tax_id=taxid(info),
        database="EXPRESSION_ATLAS",
        sequence=info["seq"].replace(
            "U", "T"
        ),  # Make sure we store the cDNA rather than RNA sequence
        regions=region_builder(info),
        rna_type=rna_type(info["rna_type"]),
        url=url(info),
        seq_version="1",
        description=info["description"],
        species=species(info),
        common_name=common_name(info),
        lineage=lineage(info),
        gene=info["gene"],
        note_data={"experiments": info["experiment"]},
    )


def find_all_taxids(directory):
    """
    Find all the taxids mentioned in all the SDRF files in the given directory.

    This will recursively glob, so should find everything within the EA cache dir.
    It also means that this one function processes all N thousand SDRF files, so
    we need a good way to catch errors here and not crash.
    """
    directory = Path(directory)
    sdrfs = []
    for p in directory.iterdir():
        if p.is_dir():
            this_sdrf = list(p.glob("*condensed-sdrf.tsv"))
            if len(this_sdrf) == 0:
                LOGGER.warning("No SDRF found in %s", p)
                continue

            sdrfs.append(this_sdrf[0])
    # sdrfs = list(directory.rglob("*condensed-sdrf.tsv"))
    sdrf_data = None
    for s in sdrfs:
        if sdrf_data is None:
            sdrf_data = sdrf.parse_condensed_sdrf(s)
        else:
            sdrf_data = pl.concat((sdrf_data, sdrf.parse_condensed_sdrf(s)))
    if sdrf_data is None:
        raise FileNotFoundError("No SDRF files found? Check the provided directory")
    organisms = (
        sdrf_data.filter(pl.col("feat_type") == "organism").select("ontology").unique()
    )
    N_not_NCBI = organisms.filter(pl.col("ontology").str.contains("NCBI").not_()).height
    ## Check if any ontology terms in the organism column are not NCBI taxa
    ## Remove them if not and give a warning
    if N_not_NCBI > 0:
        LOGGER.warning("%s taxids found that are not NCBI taxa!", N_not_NCBI)
        organisms = organisms.filter(pl.col("ontology").str.contains("NCBI"))

    try:
        taxids = (
            organisms.with_columns(
                taxid=pl.col("ontology")
                .str.split("NCBITaxon_")
                .list.last()
                .cast(pl.Int64)
            )
            .select(pl.col("taxid").unique())
            .filter(pl.col("taxid").is_not_null())
            .sort(by="taxid")
        )
    except pl.exceptions.InvalidOperationError as e:
        raise ValueError(f"Failed to extract taxids from all SDRF files: {e}")

    return taxids.get_column("taxid").to_list()
