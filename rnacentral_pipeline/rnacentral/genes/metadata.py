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

"""
Gene description and RNA type, computed without touching the database.

The evidence (accessions, Rfam and R2DT hits for the taxon's URS_taxids) is
extracted once in bulk; `compute_metadata` then scores it in polars and picks the
best description and RNA type per gene.
"""

import re
from functools import lru_cache

import networkx as nx
import obonet as obo
import polars as pl

from rnacentral_pipeline.databases.data import Database
from rnacentral_pipeline.rnacentral.precompute import utils

SO_ONTOLOGY_URL = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/heads/master/Ontology_Files/so-simple.obo"

description_scores = {
    Database.mirbase.value.descr: 100,  ## Pretty good for miRNA
    Database.wormbase.value.descr: 100,  ## Good for worm, but sometimes vague
    Database.hgnc.value.descr: 100,  ## Gold standard
    Database.gencode.value.descr: 70,  ## pretty good
    Database.ensembl.value.descr: 50,  ## Variable quality...
    Database.tair.value.descr: 100,  ## Good for arabidopsis
    Database.sgd.value.descr: 100,  ## SGD is usually pretty good
    Database.flybase.value.descr: 100,  ## Good, but sometimes vague
    Database.dictybase.value.descr: 100,  ## Very good
    Database.pombase.value.descr: 70,  ## Pombase's descriptions are a bit vague
    Database.japonicusdb.value.descr: 30,  ## Newer than the rest; a fission yeast MOD like PomBase, but lower
    Database.mirtrondb.value.descr: 20,  ## Newer; miRNA-specific but little description
    Database.psicquic.value.descr: 10,  ## Newer; interaction records, vague
    Database.ribovision.value.descr: 10,  ## Newer; structure records, vague
    Database.mgi.value.descr: 80,  ## MGI's descriptions have a lot of 'uncharacterised' in them
    Database.rgd.value.descr: 70,  ## Similar for RGD
    Database.zfin.value.descr: 50,  ## They're mostly just gene names
    Database.mirgenedb.value.descr: 90,  ## Very good for miRNA
    Database.plncdb.value.descr: 10,  ## Should prefer other plant sources if available
    Database.lncipedia.value.descr: 40,  ## Not bad but a bit vague
    Database.lncrnadb.value.descr: 90,  ## Pretty good for lncRNA
    Database.lncbook.value.descr: 40,  ## Vague, just gene names mostly
    Database.gtrnadb.value.descr: 60,  ## Just a gene name and species, but pretty good for tRNA
    Database.tmrna_website.value.descr: 90,  ## Should be the standard for tmRNA, hopefully no type conflicts
    Database.five_srrnadb.value.descr: 30,  ## Just species name + 5SrRNA, not that informative
    Database.ribocentre.value.descr: 40,  ## Can be a bit formulaic and uninformative
    Database.pdbe.value.descr: 40,  ## Can be really informative, or can be crazy long strings of nonsense
    Database.refseq.value.descr: 30,  ## Not very informative
    Database.ensembl_plants.value.descr: 15,  ## Slightly preferable to plncdb if available
    Database.ensembl_metazoa.value.descr: 10,  ## Just species + gene name, not very informative
    Database.ensembl_protists.value.descr: 20,  ## Slicghtly more informative than metazoa
    Database.ensembl_fungi.value.descr: 40,  ## Needs to be lower than SGD andPomBase, but is still quite informative
    Database.genecards.value.descr: 40,  ## All human, all homo sapiens + gene ID, so not that informative
    Database.malacards.value.descr: 40,  ## As with genecards
    Database.intact.value.descr: 40,  ## Not as good as mir-specifics, and not quite as good as tRNA specifics, but still pretty good
    Database.expression_atlas.value.descr: 10,  ## Only species + gene ID
    Database.rfam.value.descr: 65,  ## Place these slightly above others, but not by much. Now includes cm description too.
    Database.tarbase.value.descr: 10,  ## Species + gene name
    Database.lncbase.value.descr: 10,  ## Species + gene name
    Database.snodb.value.descr: 50,  ## Slightly below Rfam, but better than others
    Database.snorna_database.value.descr: 40,  ## Slightly worse than snodb
    Database.pirbase.value.descr: 0,  ## piRNAs will be getting filtered, so this should not contribute
    Database.modomics.value.descr: 10,  ## Species + type
    Database.vega.value.descr: 20,  ## Long and potentially not that valuable
    Database.srpdb.value.descr: 30,  ## Quite informative, but should be well below MODs
    Database.snopy.value.descr: 10,  ## Species + type, not that informative
    Database.crw.value.descr: 10,  ## Very short, not very informative
    Database.silva.value.descr: 5,  ## Not very informative
    Database.greengenes.value.descr: 5,  ## Not very informative
    Database.rdp.value.descr: 5,  ## Species + type
    Database.ena.value.descr: 2,  ## Often vague and uninformative
    Database.zwd.value.descr: 2,  ## Not very uninformative
    Database.noncode.value.descr: 5,  ## Species + ncRNA - not very useful
    Database.evlncrnas.value.descr: 1,  ## Very vague
    Database.mgnify.value.descr: 1,  ## All metagemonic, not super useful
    "R2DT": 30,  ## Can be very specific, but only has the model name, so not that informative
}

## Database type scoring should be slightly different
type_scores = {
    Database.five_srrnadb.value.descr: 100,
    Database.flybase.value.descr: 100,
    Database.gtrnadb.value.descr: 60,
    Database.lncbase.value.descr: 50,
    Database.lncipedia.value.descr: 40,
    Database.mirbase.value.descr: 100,
    Database.mirgenedb.value.descr: 90,
    Database.pirbase.value.descr: 50,
    Database.pombase.value.descr: 70,
    Database.snodb.value.descr: 50,
    Database.snorna_database.value.descr: 40,
    Database.tarbase.value.descr: 100,
    Database.zwd.value.descr: 2,
    Database.tmrna_website.value.descr: 90,
    Database.pdbe.value.descr: 40,
    Database.gencode.value.descr: 70,
    Database.genecards.value.descr: -20,  ## These should be ignored
    Database.malacards.value.descr: -20,
    Database.rfam.value.descr: 100,
    "R2DT": 100,  ## Can be very specific
}


@lru_cache()
def get_so_graph():
    # Lazy: every `rnac` command imports this module (cli/genes -> classify ->
    # here), so an import-time fetch broke commands that never use the ontology.
    return obo.read_obo(SO_ONTOLOGY_URL)


# Generic "ncRNA" SO term. Any SO id we can't find in the loaded ontology
# (renamed/removed across SO releases, or never present) gets short-circuited
# to this so we degrade to crappy-but-valid typing instead of erroring out.
GENERIC_NCRNA = "SO:0000655"


def normalize_so_type(so_type):
    """Map an SO id to itself if it exists in the ontology, else GENERIC_NCRNA.

    None and any id that is not a node in so_graph collapse to the generic
    ncRNA term, which is guaranteed to be present.
    """
    if so_type is not None and so_type in get_so_graph():
        return so_type
    return GENERIC_NCRNA


def calculate_type_specificity(so_type):
    # Unknown/missing SO ids are short-circuited to the generic ncRNA term
    # rather than blowing up (NodeNotFound) or being heavily penalised.
    so_type = normalize_so_type(so_type)
    try:
        return nx.shortest_path_length(get_so_graph(), so_type, "SO:0000673")
    except (nx.NetworkXNoPath, nx.NodeNotFound):
        return -1000


def remove_species_mention(description, species_name=None, common_name=None):
    """
    Remove species mentions from descriptions to make a short description.

    Args:
        description: The text to clean
        species_name: Scientific name (e.g., "Homo sapiens") - optional
        common_name: Common name (e.g., "human") - optional
    """

    # If specific species/common names are provided, remove them
    if species_name:
        # Remove scientific name (case insensitive)
        description = re.sub(
            rf"\b{re.escape(species_name)}\b", "", description, flags=re.IGNORECASE
        )

    if common_name:
        # Remove common name in parentheses
        description = re.sub(
            rf"\s*\(\s*{re.escape(common_name)}\s*\)",
            "",
            description,
            flags=re.IGNORECASE,
        )
        # Also remove standalone common name
        description = re.sub(
            rf"\b{re.escape(common_name)}\b", "", description, flags=re.IGNORECASE
        )

    # General pattern to catch common formats if no specific names provided
    if not species_name and not common_name:
        # Remove pattern like "Genus species (common name)"
        description = re.sub(r"\b[A-Z][a-z]+ [a-z]+ \([^)]+\)", "", description)

    # Clean up extra whitespace
    description = re.sub(r"\s+", " ", description).strip()

    return description


EVIDENCE_SCHEMA = {
    "urs_taxid": pl.Utf8,
    "database": pl.Utf8,
    "description": pl.Utf8,
    "rna_type": pl.Utf8,
    "cm_overlap": pl.Float64,
}

METADATA_SCHEMA = {
    "internal_name": pl.Utf8,
    "description": pl.Utf8,
    "rna_type": pl.Utf8,
    "short_description": pl.Utf8,
}


def _lookup(values, func):
    """Apply func once per distinct value, then map back (entropy and graph walks are slow)."""
    mapping = {v: func(v) for v in values.drop_nulls().unique().to_list()}
    return mapping


def compute_metadata(gene_members, evidence):
    """
    gene_members - frame with internal_name and region_name (one row per member)
    evidence     - frame with EVIDENCE_SCHEMA

    Returns (metadata, missing): metadata has METADATA_SCHEMA, one row per gene that
    has any evidence; missing lists the genes that have none. Ties are broken on the
    description / rna_type text so the result does not depend on row order.
    """
    pairs = gene_members.select(
        "internal_name",
        pl.col("region_name").str.split("@").list.first().alias("urs_taxid"),
    ).unique()
    ev = evidence.cast(EVIDENCE_SCHEMA).with_columns(
        pl.col("cm_overlap").fill_null(0.0), pl.col("description").fill_null("")
    )

    entropy = _lookup(ev["description"], utils.entropy)
    specificity = _lookup(ev["rna_type"], calculate_type_specificity)
    generic = normalize_so_type(None)

    ev = ev.with_columns(
        (
            pl.col("database").replace_strict(description_scores, default=0)
            + pl.col("description").replace_strict(entropy, default=0.0)
        ).alias("desc_score"),
        (
            pl.col("database").replace_strict(type_scores, default=0)
            + pl.col("rna_type").replace_strict(specificity, default=-1000)
            + pl.col("cm_overlap") * 100.0
        ).alias("type_score"),
    )
    joined = pairs.join(ev, on="urs_taxid", how="inner")

    best_desc = (
        joined.sort(
            ["internal_name", "desc_score", "description"],
            descending=[False, True, False],
        )
        .unique(subset="internal_name", keep="first", maintain_order=True)
        .select("internal_name", "description")
    )
    best_type = (
        joined.sort(
            ["internal_name", "type_score", "rna_type"],
            descending=[False, True, False],
            nulls_last=True,
        )
        .unique(subset="internal_name", keep="first", maintain_order=True)
        .select("internal_name", "rna_type")
    )
    metadata = best_desc.join(best_type, on="internal_name")
    ## A missing or unknown SO id collapses to the generic ncRNA term
    normalized = _lookup(metadata["rna_type"], normalize_so_type)
    shorts = _lookup(metadata["description"], remove_species_mention)
    metadata = metadata.with_columns(
        pl.col("rna_type").replace_strict(normalized, default=generic),
        pl.col("description")
        .replace_strict(shorts, default="")
        .alias("short_description"),
    ).sort("internal_name")

    missing = (
        pairs.select("internal_name")
        .unique()
        .join(metadata.select("internal_name"), on="internal_name", how="anti")
        .sort("internal_name")
    )
    return metadata.select(list(METADATA_SCHEMA)), missing
