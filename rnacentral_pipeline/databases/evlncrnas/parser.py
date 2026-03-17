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
from rnacentral_pipeline.databases.data import Entry, Exon, SequenceRegion
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.rnacentral import lookup

from . import helpers

tqdm.pandas()

QUERY = """
select
    pre.id as id,
    pre.rna_type,
    COALESCE(rna.seq_short, rna.seq_long) as sequence,
    pre.description

from rnc_rna_precomputed pre
join rna on rna.upi = pre.upi
where
    pre.id in %s
"""

base_url = furl(
    "https://www.sdklab-biophysics-dzu.net/EVLncRNAs2/index.php/Home/Browsc/rna.html"
)

ensembl_rest_url = furl("https://rest.ensembl.org/sequence/id")

entrez_base_url = furl("https://eutils.ncbi.nlm.nih.gov/entrez/eutils")

chain_normalisation = {
    "minus": "-",
    "plus": "+",
}

type_normalisation = {
    "lncRNA": "SO:0001877",
    "lincRNA": "SO:0001463",
    "antisense": "SO:0001904",
}


def handled_phylogeny(species: str) -> int:
    try:
        return phy.taxid(species)
    except phy.FailedTaxonId:
        return None


def condense_publications(record):
    references = []
    seen = set()
    for value in record:
        if pd.isna(value):
            continue
        try:
            pmid = int(value)
        except (TypeError, ValueError):
            continue
        if pmid in seen:
            continue
        seen.add(pmid)
        references.append(pubs.reference(pmid))
    return references


def resolve_sheet(db_dir: Path, basename: str) -> Path:
    for suffix in (".xls", ".tsv"):
        candidate = db_dir.joinpath(f"{basename}{suffix}")
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find {basename}.xls or {basename}.tsv in {db_dir}")


def load_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    return pd.read_excel(path)


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
    e_accessions = no_accessions[no_accessions["Ensembl"].notna()].copy()
    print("ensembl subset done")
    no_accessions = no_accessions[no_accessions["Ensembl"].isna()].copy()
    ncbi_accessions = input_frame[input_frame["NCBI accession"].notna()].copy()
    print("NCBI subset done")
    return (no_accessions, e_accessions, ncbi_accessions)


def get_ncbi_accessions(
    accession_frame_in: pd.DataFrame,
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
                fetch_url.args["rettype"] = "gb"
                fetch_url.args["retmode"] = "text"

                fetch_res = requests.get(fetch_url.url)
                if fetch_res.ok:
                    sequence_data = SeqIO.parse(
                        io.StringIO(fetch_res.text),
                        "gb",
                    )
                    for record in sequence_data:
                        sequences.append(str(record.seq).replace("U", "T"))
        return sequences

    accession_frame = accession_frame_in.copy()
    accession_frame["sequence"] = accession_frame["NCBI accession"].apply(
        download_and_get_sequence
    )
    accession_frame = accession_frame.explode("sequence")
    missing_frame = accession_frame[accession_frame["sequence"].isna()]

    return (accession_frame.dropna(subset="sequence"), missing_frame)


def get_ensembl_accessions(
    ensembl_frame_in: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def pull_ensembl_data(e_id: str):
        id_url = ensembl_rest_url / e_id

        data = requests.get(id_url.url, headers={"Content-Type": "text/x-fasta"})
        if not data.ok:
            return (None, None, None, None, None, None)
        sequence_data = SeqIO.read(
            io.StringIO(data.text),
            "fasta",
        )
        details = sequence_data.description.split(":")
        assembly = details[1]
        chromosome = details[2]
        region_start = details[3]
        region_stop = details[4]
        strand = details[5]

        sequence = str(sequence_data.seq).replace("U", "T")
        return (sequence, assembly, chromosome, region_start, region_stop, strand)

    ensembl_frame = ensembl_frame_in.copy()
    ensembl_frame[
        [
            "sequence",
            "assembly_id",
            "chromosome",
            "region_start",
            "region_stop",
            "chain",
        ]
    ] = ensembl_frame.apply(
        lambda row: pull_ensembl_data(row["Ensembl"]),
        axis="columns",
        result_type="expand",
    )
    missing_frame = ensembl_frame[ensembl_frame["sequence"].isna()]
    return (ensembl_frame.dropna(subset="sequence"), missing_frame)


def get_db_matches(match_frame_in: pd.DataFrame, db_dump: Path) -> pd.DataFrame:
    def lookup_names(row):
        names = [str(row["external_id"]).strip()]
        aliases = row.get("Aliases")
        if pd.notna(aliases):
            names.extend(a.strip() for a in str(aliases).split(",") if a.strip())
        return names

    match_frame = match_frame_in.copy()
    match_frame["taxid"] = match_frame["taxid"].astype(int)
    match_frame["lookup_name"] = match_frame.apply(lookup_names, axis="columns")
    match_frame = (
        match_frame.explode("lookup_name")
        .replace(to_replace=["None"], value=np.nan)
        .dropna(subset="lookup_name")
    )
    match_frame["is_exact_match"] = (
        match_frame["lookup_name"] == match_frame["external_id"]
    )

    rnc_data = pd.read_csv(db_dump, names=["urs", "taxid", "lookup_name"], header=0)
    rnc_data["lookup_name"] = rnc_data["lookup_name"].apply(lambda x: str(x).split("|"))
    rnc_data = (
        rnc_data.explode("lookup_name")
        .replace(to_replace=["", None], value=np.nan)
        .dropna(subset="lookup_name")
    )

    matches = match_frame.merge(
        rnc_data,
        left_on=["lookup_name", "taxid"],
        right_on=["lookup_name", "taxid"],
        how="inner",
    )
    matches.sort_values(["ID", "is_exact_match"], ascending=[True, False], inplace=True)

    return matches


def load_function_data(function_info: Path) -> pd.DataFrame:
    function_df = load_table(function_info)
    return (
        function_df.groupby("ID", sort=False)["PMID"]
        .apply(condense_publications)
        .reset_index(name="publications")
    )


def parse(db_dir: Path, db_dumps: tuple[Path], db_url: str) -> None:
    """
    Parse and join the two EVLncRNAs3 workbooks and build RNAcentral entries.
    """
    lncRNA = resolve_sheet(db_dir, "lncRNA_information")
    function_info = resolve_sheet(db_dir, "function_information")

    lncRNA_df = load_table(lncRNA)
    function_df = load_function_data(function_info)
    lncRNA_df.rename(
        columns={
            "LncRNA name": "external_id",
            "Alias": "Aliases",
        },
        inplace=True,
    )

    print("Loaded EVLncRNAs3 sheets...")

    lncRNA_df["taxid"] = lncRNA_df["Species"].apply(handled_phylogeny)
    lncRNA_df = lncRNA_df.dropna(subset=["taxid"]).copy()
    lncRNA_df["taxid"] = lncRNA_df["taxid"].astype(int)

    ## Split the data on the presence of accessions for either NCBI or Ensembl
    no_accession_frame, ensembl_frame, ncbi_frame = split(lncRNA_df)

    ## These two look up directly from the source, so should be quick ish
    ensembl_frame, missing_ensembl_frame = get_ensembl_accessions(ensembl_frame)
    print(f"Got all available ensembl accessions ({len(ensembl_frame)})")

    ncbi_frame, missing_ncbi_frame = get_ncbi_accessions(ncbi_frame)
    print(f"Got all available NCBI accessions ({len(ncbi_frame)})")

    ## Stack the frames with missing accessions together to search RNAcentral
    no_accession_frame = pd.concat(
        [no_accession_frame, missing_ensembl_frame, missing_ncbi_frame]
    )  #

    ## Match with RNAcentral based on the gene name
    ## This is optionally chunked to save memory - 
    ## split the lookup file and provide a list on the commandline
    matched_chunks = [get_db_matches(no_accession_frame, dump_chunk) for dump_chunk in db_dumps]
    matched_frame = pd.concat(matched_chunks, ignore_index=True)
    matched_frame.drop_duplicates(subset="ID", inplace=True)
    matched_frame["urs_taxid"] = (
        matched_frame["urs"] + "_" + matched_frame["taxid"].astype(str)
    )

    ## Look up the rest of the data for the hits
    mapping = lookup.as_mapping(db_url, matched_frame["urs_taxid"].values, QUERY)
    for idx, value in enumerate(mapping.values()):
        value["sequence"] = value["sequence"].replace("U", "T")

    ## Copy matching sequence data
    matched_frame["sequence"] = matched_frame["urs_taxid"].apply(
        lambda x: mapping[x]["sequence"]
    )

    ## Build frame with all hits & accessions and add aggregated publication data
    full_frame = pd.concat([matched_frame, ensembl_frame, ncbi_frame], ignore_index=True)
    full_frame.drop_duplicates(subset="ID", inplace=True)
    full_frame = full_frame.merge(function_df, how="left", on="ID")
    full_frame["publications"] = full_frame["publications"].apply(
        lambda refs: refs if isinstance(refs, list) else []
    )

    ## Tidy up and apply some normalisations
    full_frame["Chain"] = full_frame["Chain"].apply(
        lambda x: chain_normalisation.get(str(x).lower(), None) if pd.notna(x) else None
    )
    full_frame["so_type"] = full_frame["Class"].apply(
        lambda x: type_normalisation.get(x, "SO:0000655")
    )

    full_frame.replace({np.nan: None}, inplace=True)

    ## yield entry objects for each row in the frame, these get written directly.
    for _, row in full_frame.iterrows():
        yield helpers.as_entry(row)
