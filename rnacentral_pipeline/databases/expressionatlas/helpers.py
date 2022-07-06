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

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs


def accession(urs_taxid):
    return "EXPRESSIONATLAS:" + urs_taxid


def primary_id(urs_taxid):
    return "EXPRESSIONATLAS:" + urs_taxid


def taxid(urs_taxid):
    _, taxid = urs_taxid.split("_")
    return int(taxid)


def species(urs_taxid):
    return phy.species(taxid(urs_taxid))


def lineage(urs_taxid):
    return phy.lineage(taxid(urs_taxid))


def common_name(urs_taxid):
    return phy.common_name(taxid(urs_taxid))


def url(ensembl_id):
    return "https://www.ebi.ac.uk/gxa/genes/" + ensembl_id


def references(interactions):
    refs = set()
    for interaction in interactions:
        refs.update(interaction.publications)
    refs.add(pubs.reference(24234451))
    return list(refs)



def as_entry(urs_taxid, info, gene):
    return data.Entry(
        primary_id=primary_id(urs_taxid),
        accession=accession(urs_taxid),
        ncbi_tax_id=taxid(urs_taxid),
        database="EXPRESSIONATLAS",
        sequence=info["sequence"],
        regions=[],
        rna_type=info["rna_type"],
        url=url(gene),
        seq_version="1",
        description=info["description"],
        species=species(urs_taxid),
        common_name=common_name(urs_taxid),
        lineage=lineage(urs_taxid),
    )
