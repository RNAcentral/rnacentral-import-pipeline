# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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

import operator as op
from functools import partial
from operator import is_not
from pathlib import Path

import numpy as np
import pandas as pd
from furl import furl

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy

from . import lookup

base_url = furl(
    "https://www.sdklab-biophysics-dzu.net/EVLncRNAs2/index.php/Home/Browsc/rna.html"
)


def handled_phylogeny(species):
    try:
        return phy.taxid(species)
    except phy.FailedTaxonId:
        return None


def as_entry(record):
    """
    Generate an Entry to import based off the database, exons and raw record.
    """

    if len(record["exons"]) == 0:
        print(
            f"Something is wrong with sequence {record['index']}, the exon list is empty"
        )
        chromosome = str(record["chromosome"])
        if record["Chain"] is not None:
            strand = (
                -1 if record["Chain"] == "minus" else 1
            )  ## This may not be right... think mito RNAs?
        else:
            strand = None
    else:
        chromosome = record["exons"][0]["chromosome"]
        strand = record["exons"][0]["strand"]

    region = data.SequenceRegion(
        chromosome=chromosome,
        strand=strand,
        exons=[
            data.Exon(start=e["exon_start"], stop=e["exon_stop"])
            for e in record["exons"]
        ],
        assembly_id=record["Assembly"],
        coordinate_system="1-start, fully-closed",
    )
    try:
        entry = data.Entry(
            primary_id=record["index"],
            accession=record["index"],
            ncbi_tax_id=record["taxid"],
            database="EVLNCRNA",
            sequence=record["sequence"],
            regions=[region],
            gene=record.get("Name", ""),
            gene_synonyms=record.get("Aliases", []),
            rna_type=record.get("so_term", record["Class"]),
            url=base_url.set({"id": record["index"]}).url,
            seq_version=record.get("version", "1"),
            optional_id=record.get("optional_id", ""),
            description=record.get("description", ""),
            # note_data=note_data(record),
            # xref_data=xrefs(record),
            # related_sequences=related_sequences(record),
            species=record.get("Species", ""),
            lineage=phy.lineage(record["taxid"]),
            common_name=phy.common_name(record["taxid"]),
            chromosome=str(record["chromosome"]),
            # secondary_structure=secondary_structure(record),
            # references=references(record),
            organelle=record.get("localization", None),
            product=record.get("product", None),
            # location_start=record.get("Start site", None),
            # location_end=record.get("End site", None),
            # anticodon=anticodon(record),
            # gene=gene(record),
            # gene_synonyms=gene_synonyms(record),
            # locus_tag=locus_tag(record),
            # features=features(record),
        )
    except:
        print(
            f"Unexpected RNA type: {record['so_term']} for record with ID {record['index']}"
        )
        entry = None

    return entry


def parse(db_dir: Path, db_url: str):
    """
    Parses the 3 excel sheets using pandas and joins them into one massive table
    which is then parsed to produce entries
    """
    lncRNA = db_dir.joinpath("lncRNA.xlsx")
    interaction = db_dir.joinpath("interaction2.xlsx")
    disease = db_dir.joinpath("disease2.xlsx")

    assert lncRNA.exists() and interaction.exists() and disease.exists()

    lncRNA_df = pd.read_excel(lncRNA, index_col=0)
    interaction_df = pd.read_excel(interaction, index_col=0)
    disease_df = pd.read_excel(disease, index_col=0)

    print("Loaded 3 sheets...")

    ## Joins on index if no on= is specified. The index is the ID column because
    ## we set index_col=0 during load, so this is same as saying on='ID'
    ## The benefit of doing it this way is that we can join multiple dfs on the
    ## ID because it is index
    complete_df = lncRNA_df.join(
        [
            disease_df.drop(
                columns=["Name", "Species", "Species category", "exosome", "structure"]
            ),
            interaction_df.drop(columns=["Name", "Species", "Species category"]),
        ],
        how="outer",
        sort=True,
    )

    print("Sheet join complete...")

    ## Fix the way NA is read as NaN, cast to string
    ## Empty strings are dealt with in the lookup.mapping function
    complete_df.replace(np.nan, "", inplace=True)

    ## Expand the list of aliases by preprending the Name of the RNA, then splitting on commas
    complete_df["Aliases"] = (
        complete_df["Name"] + ", " + complete_df["Aliases"].astype(str)
    ).str.split(",")

    complete_df["taxid"] = complete_df["Species"].apply(handled_phylogeny)

    ## Lookup the aliases to get what data we hold about them. This is
    ## ext_id -> RNAc information about it
    # mapping = lookup.mapping(db_url, {ID:aliases for ID, aliases in zip(complete_df.index.values, complete_df["Aliases"].values)})
    print("Getting ID lookup from database...")
    mapping = lookup.mapping(db_url, complete_df)

    # complete_df = complete_df.join(mapping, how='outer', sort=True)
    complete_df.replace(np.nan, "", inplace=True)
    # complete_df = complete_df[complete_df['urs_taxid'] != ""]

    print(complete_df.groupby("Species").count())
    ## Add a column with taxid looked up from the species

    complete_df = complete_df.reset_index()

    entries = complete_df.apply(as_entry, axis=1)

    return list(filter(partial(is_not, None), entries))
