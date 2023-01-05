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
import subprocess as sp
import zipfile
from functools import partial
from operator import is_not
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from Bio import SeqIO
from furl import furl
from tqdm import tqdm

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy

from . import lookup

tqdm.pandas()

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


def split(db_dir: Path, output_loc: Path):
    """
    Split the main dataset based on presence of an NCBI accession
    """
    lncRNA = db_dir.joinpath("lncRNA.xlsx")

    no_acc = output_loc.joinpath("no_accessions.jsonl")
    acc = output_loc.joinpath("with_accessions.jsonl")

    assert lncRNA.exists()
    lncRNA_df = pd.read_excel(lncRNA)

    ## We might as well lookup the taxid while we're here
    lncRNA_df["taxid"] = lncRNA_df["Species"].apply(
        handled_phylogeny
    )  # .dropna().astype(int)

    no_accessions = lncRNA_df[lncRNA_df["NCBI accession"].isna()]
    accessions = lncRNA_df[lncRNA_df["NCBI accession"].notna()]
    print(len(no_accessions), len(accessions))

    with open(no_acc, "w") as no_acc_output:
        no_acc_output.write(
            f"{no_accessions[['ID', 'Name', 'taxid']].to_json(orient='records', lines=True)}"
        )
    with open(acc, "w") as acc_output:
        acc_output.write(
            f"{accessions[['ID', 'Name', 'taxid', 'NCBI accession']].to_json(orient='records', lines=True)}"
        )


def get_accessions(accession_file: Path, output: Path):
    """
    For each entry having at least one NCBI accession, build the commandline for the NCBI
    datasets tool and download the data. Then extract the RNA sequence and build a new
    dataframe with the necessary information.
    """

    def download_and_get_sequence(accessions):
        temp_filename = NamedTemporaryFile(delete=False).name
        accession_list = [x.strip() for x in accessions.split(",")]
        command = base_command.format(
            accession_list[0], temp_filename, ",".join(accession_list)
        ).split()
        sp.Popen(command, stdin=None, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        try:
            downloaded_data = zipfile.ZipFile(temp_filename, "r")
        except zipfile.BadZipFile:
            if len(accession_list) == 1:
                print(
                    f"No sequences found for accession {accession_list}, returning empty list"
                )
                return []
            else:
                """
                try to recover by moving to the next accession
                """
                print(accession_list)
                accession_list.pop(0)
                print(accession_list)
                return download_and_get_sequence(",".join(accession_list))

        sequence_data = SeqIO.parse(
            io.TextIOWrapper(
                downloaded_data.open("ncbi_dataset/data/rna.fna", "r"), encoding="utf-8"
            ),
            "fasta",
        )
        sequences = []
        for record in sequence_data:
            sequences.append("".join(record.seq).replace("U", "T"))
        return sequences

    assert accession_file.exists()
    lnc_data = pd.read_json(accession_file, lines=True)
    print(lnc_data)

    base_command = "bin/datasets download gene accession {0} --filename {1} --include rna --fasta-filter {2} --no-progressbar"
    lnc_data["sequences"] = lnc_data["NCBI accession"].progress_apply(
        download_and_get_sequence
    )
    lnc_data = lnc_data.explode("sequences").dropna()

    with open(output, "w") as out_file:
        out_file.write(f"{lnc_data.to_json(orient='records', lines=True)}")


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
    complete_df["taxid"] = complete_df["Species"].apply(handled_phylogeny).astype(int)

    print(complete_df)

    exit()
