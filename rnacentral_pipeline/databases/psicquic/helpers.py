# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import typing as ty

from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.data import Entry, Interaction


def taxid(urs_taxid):
    _, taxid = urs_taxid.split("_")
    return int(taxid)


def species(interaction):
    return phy.species(taxid(interaction))


def lineage(interaction):
    return phy.lineage(taxid(interaction))


def common_name(interaction):
    return phy.common_name(taxid(interaction))


def as_entry(urs_taxid: str, interactions: ty.List[Interaction], info: ty.Dict[str, str]):
    database = "PSICQUIC"

    return Entry(
        primary_id=f"{database}:{urs_taxid}",
        accession=f"{database}:{urs_taxid}",
        ncbi_tax_id=taxid(urs_taxid),
        database=database,
        sequence=info["sequence"],
        regions=[],
        rna_type=info["rna_type"],
        url=url(urs_taxid),
        seq_version="1",
        description=info["description"],
        references=references(interactions),
        species=species(urs_taxid),
        common_name=common_name(urs_taxid),
        lineage=lineage(urs_taxid),
        interactions=interactions,
    )
