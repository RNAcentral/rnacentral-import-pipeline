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
from rnacentral_pipeline.databases.helpers import publications as pub
from rnacentral_pipeline.databases.data import Entry, Interaction, IdReference


def taxid(urs_taxid: str):
    _, taxid = urs_taxid.split("_")
    return int(taxid)


def species(urs_taxid: str):
    return phy.species(taxid(urs_taxid))


def lineage(urs_taxid: str):
    return phy.lineage(taxid(urs_taxid))


def common_name(urs_taxid: str):
    return phy.common_name(taxid(urs_taxid))


def references(interactions: ty.List[Interaction]) -> ty.List[IdReference]:
    ids = [
        pub.reference('PMID:23671334')
    ]
    for interaction in interactions:
        ids.extend(interaction.publications)
    return ids


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
        url="http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml",
        seq_version="1",
        description=info["description"],
        references=references(interactions),
        species=species(urs_taxid),
        common_name=common_name(urs_taxid),
        lineage=lineage(urs_taxid),
        interactions=interactions,
    )
