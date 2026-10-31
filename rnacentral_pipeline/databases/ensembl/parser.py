# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import re
import typing as ty
from pathlib import Path

import attr

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ensembl import gff
from rnacentral_pipeline.databases.ensembl.vertebrates import parser as vertebrates

# Ensembl now serves every organism (all former divisions plus bacteria) in one
# uniform EMBL/GFF3 format, so a single parser handles them all and every entry
# is imported as the ENSEMBL database. A couple of per-organism corrections from
# the old per-division parsers are preserved below.

URL = "https://www.ensembl.org/feature-explorer/{accession}/transcript:{transcript}"


def correct_protist_rna_type(entry: Entry) -> Entry:
    """
    Leishmania sno/snRNA genes are mislabelled in the source; fix them by gene
    name and description. A no-op for every other organism.
    """
    if re.match(r"^LMJF_\d+_snoRNA_?\d+$", entry.gene):
        return attr.evolve(entry, rna_type="snoRNA")
    if re.match(r"^LMJF_\d+_snRNA_\d+$", entry.gene):
        return attr.evolve(entry, rna_type="snRNA")
    if entry.rna_type == "snRNA" and "snoRNA" in entry.description:
        return attr.evolve(entry, rna_type="snoRNA")
    return entry


def as_tair_entry(entry: Entry) -> Entry:
    database = "TAIR"
    xrefs = dict(entry.xref_data)
    xrefs.pop(database, None)
    return attr.evolve(
        entry,
        accession="%s:%s" % (database, entry.primary_id),
        database=database,
        xref_data=xrefs,
    )


def tair_entries(entry: Entry) -> ty.Iterable[Entry]:
    """
    Arabidopsis thaliana ncRNAs are additionally attributed to TAIR, which
    RNAcentral tracks as its own database. This mirrors the old Ensembl Plants
    behaviour (the species JSON labels the provider only as "community").
    """
    if entry.ncbi_tax_id == 3702 and entry.primary_id.startswith("AT"):
        yield as_tair_entry(entry)


def parse(
    raw: ty.IO, gff_file: Path, family_file=None, excluded_file=None
) -> ty.Iterable[Entry]:
    accession = gff.get_assembly_accession(gff_file)
    entries = vertebrates.parse(
        raw, gff_file, family_file=family_file, excluded_file=excluded_file
    )
    for entry in entries:
        entry = correct_protist_rna_type(entry)
        entry = attr.evolve(
            entry, url=URL.format(accession=accession, transcript=entry.primary_id)
        )
        yield entry
        yield from tair_entries(entry)
