# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import csv
import itertools as it
import operator as op
import typing as ty
from pathlib import Path

from Bio import SeqIO

from rnacentral_pipeline.databases.data import Entry, Interaction
from rnacentral_pipeline.databases.helpers import phylogeny as phy


def cleaned(raw: ty.Dict[str, str]) -> ty.Dict[str, None | str]:
    """Replace 'NA' values with None in a dict.

    >>> cleaned({"a": "1", "b": "NA"})
    {"a": "1", "b": None}
    >>> cleaned({"a": "1", "b": "bob"})
    {"a": "1", "b": "bob"}
    """

    entry = {}
    for key, value in raw.items():
        val = value
        if val == "NA":
            val = None
        entry[key] = val
    return entry


def build_interactions(raw: ty.List[ty.Dict[str, str]]) -> ty.List[Interaction]:
    data = [cleaned(r) for r in raw]
    raise ValueError("Implement me")


def parse(handle: ty.IO, fasta: Path) -> ty.Iterable[Entry]:
    """Parse the given TSV file handle and FASTA file path and produce an
    iterable of Entry objects. This assumes that the file is from Tarbase v9
    and is sorted by the mirna_name column. It will produce a single Entry
    object per mirna_name with all interactions for that miRNA.
    """

    indexed = SeqIO.index(fasta, "fasta")

    reader = csv.DictReader(handle)
    grouped = it.groupby(reader, op.itemgetter("mirna_name"))
    for accession, raw_interactions in grouped:
        raw_interactions = list(raw_interactions)
        accession = raw_interactions[0]["mirna_name"]
        species = raw_interactions[0]["species"]

        sequence = indexed[accession]
        if not sequence:
            raise ValueError(f"Failed to find sequence for {accession}")
        sequence = str(sequence.seq)

        yield Entry(
            primary_id=f"MIRTARBASE:{accession}",
            accession=accession,
            ncbi_tax_id=phy.taxid(species),
            database="MIRTARBASE",
            sequence=sequence,
            regions=[],
            rna_type="SO:0000276",
            url=f"https://dianalab.e-ce.uth.gr/tarbasev9/interactions?gene=&mirna={accession}",
            seq_version="1",
            interactions=build_interactions(raw_interactions),
            description=f"{species} {accession} miRNA",
        )
