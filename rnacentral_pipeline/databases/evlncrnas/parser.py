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

import io
import operator as op
import re
from functools import partial
from operator import is_not
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from bs4 import BeautifulSoup
from furl import furl
from tqdm import tqdm

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy

from . import lookup

tqdm.pandas()

base_url = furl(
    "https://www.sdklab-biophysics-dzu.net/EVLncRNAs2/index.php/Home/Browsc/rna.html"
)

ensembl_rest_url = furl("https://rest.ensembl.org/sequence/id")

entrez_base_url = furl("https://eutils.ncbi.nlm.nih.gov/entrez/eutils")


def handled_phylogeny(species: str) -> int:
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


def split(input_frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Split the main dataset based on presence of an NCBI accession, or ensembl accession

    This will return a tuple of 3 frames, which can be dispatched to the 3 handlers
    """
    print("Splitting based on presence of accessions...")
    no_accessions = input_frame[input_frame["NCBI accession"].isna()].dropna(
        subset="taxid"
    )
    print("NCBI missing done")
    e_accessions = no_accessions[no_accessions["Ensembl"].notna()]
    print("ensembl subset done")
    ncbi_accessions = input_frame[input_frame["NCBI accession"].notna()]
    print("NCBI subset done")
    return (no_accessions, e_accessions, ncbi_accessions)


def get_ncbi_accessions(
    accession_frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each entry having at least one NCBI accession, build the commandline for the NCBI
    datasets tool and download the data. Then extract the RNA sequence and build a new
    dataframe with the necessary information.
    """

    def download_and_get_sequence(accessions):
        accession_list = [x.strip() for x in accessions.split(",")]
        sequences = []
        search_url = entrez_base_url / "esearch.fcgi"
        search_url.args["db"] = "nuccore"
        search_url.args["term"] = " OR ".join(accession_list)
        search_url.args["usehistory"] = "y"

        search_result = requests.get(search_url.url)
        if search_result.ok:
            search_res_data = BeautifulSoup(search_result.text, features="xml")
            num_hits = search_res_data.find("Count")
            if num_hits and int(num_hits.text) > 0:
                query_key = search_res_data.find("QueryKey").text
                webenv = search_res_data.find("WebEnv")

                fetch_url = entrez_base_url / "efetch.fcgi"
                fetch_url.args["db"] = "nuccore"
                fetch_url.args["term"] = " OR ".join(accession_list)
                fetch_url.args["query_key"] = query_key
                fetch_url.args["WebEnv"] = webenv
                fetch_url.args["rettype"] = "fasta"
                fetch_url.args["retmode"] = "text"

                fetch_res = requests.get(fetch_url.url)
                if fetch_res.ok:
                    sequence_data = SeqIO.parse(
                        io.StringIO(fetch_res.text),
                        "fasta",
                    )
                    for record in sequence_data:
                        sequences.append(str(record.seq).replace("U", "T"))
        return sequences

    accession_frame["sequence"] = accession_frame["NCBI accession"].progress_apply(
        download_and_get_sequence
    )
    accession_frame = accession_frame.explode("sequence")
    missing_frame = accession_frame[accession_frame["sequence"].isna()]

    return (accession_frame.dropna(subset="sequence"), missing_frame)


def get_ensembl_accessions(
    ensembl_frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def pull_ensembl_data(e_id: str):
        id_url = ensembl_rest_url / e_id

        data = requests.get(id_url.url, headers={"Content-Type": "text/x-fasta"})
        if not data.ok:
            return None
        sequence_data = SeqIO.read(
            io.StringIO(data.text),
            "fasta",
        )
        sequence = str(sequence_data.seq).replace("U", "T")
        return sequence

    ensembl_frame["sequence"] = ensembl_frame["Ensembl"].progress_apply(
        pull_ensembl_data
    )
    missing_frame = ensembl_frame[ensembl_frame["sequence"].isna()]
    return (ensembl_frame.dropna(subset="sequence"), missing_frame)


def get_db_matches(match_frame: pd.DataFrame, db_dump: Path) -> pd.DataFrame:
    def split_clean_aliases(al):
        if al:
            return [a.strip() for a in str(al).split(",")]
        return np.nan

    print(len(match_frame))
    match_frame["taxid"] = match_frame["taxid"].astype(int)
    match_frame["external_id"] = match_frame[["Name", "Aliases"]].apply(
        lambda x: ",".join(x.values.astype(str)), axis=1
    )
    match_frame["external_id"] = match_frame["external_id"].apply(split_clean_aliases)
    match_frame = (
        match_frame.explode("external_id")
        .replace(to_replace=["None"], value=np.nan)
        .dropna(subset="external_id")
    )
    # lnc_data = lnc_data.set_index("external_id")

    rnc_data = pd.read_csv(db_dump)
    rnc_data["external_id"] = rnc_data["external_id"].apply(lambda x: str(x).split("|"))
    rnc_data = (
        rnc_data.explode("external_id")
        .replace(to_replace=["", None], value=np.nan)
        .dropna(subset="external_id")
    )
    # rnc_data = rnc_data.set_index("external_id")

    matches = match_frame.merge(
        rnc_data,
        left_on=["external_id", "taxid"],
        right_on=["external_id", "taxid"],
        how="inner",
    ).drop_duplicates(subset="urs")

    return matches


def parse(db_dir: Path, db_dump: Path, db_url: str) -> None:
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
    complete_df = lncRNA_df
    # .join(
    #     [
    #         disease_df.drop(
    #             columns=["Name", "Species", "Species category", "exosome", "structure"]
    #         ),
    #         interaction_df.drop(columns=["Name", "Species", "Species category"]),
    #     ],
    #     how="inner",
    #     sort=True,
    # )

    print("Sheet join complete...")

    ## Fix the way NA is read as NaN, cast to string
    ## Empty strings are dealt with in the lookup.mapping function
    # complete_df.replace(np.nan, "", inplace=True)
    complete_df["taxid"] = (
        complete_df["Species"].progress_apply(handled_phylogeny).dropna().astype(int)
    )

    no_accession_frame, ensembl_frame, ncbi_frame = split(complete_df)

    print(len(ensembl_frame), len(ncbi_frame), len(no_accession_frame))
    ## These two look up directly from the source, so should be quick ish
    # ensembl_frame, missing_ensembl_frame = get_ensembl_accessions(ensembl_frame)
    print(f"Got all available ensembl accessions ({len(ensembl_frame)})")

    # ncbi_frame, missing_ncbi_frame = get_ncbi_accessions(ncbi_frame)
    # print(f"Got all available NCBI accessions ({len(ncbi_frame)})")
    # no_accession_frame = pd.concat([no_accession_frame, missing_ensembl_frame]) # missing_ncbi_frame

    no_accession_frame = get_db_matches(no_accession_frame, db_dump)
    print(len(ensembl_frame), len(no_accession_frame))  # len(ncbi_frame),

    exit()
