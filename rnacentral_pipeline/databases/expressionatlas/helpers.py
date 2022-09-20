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

from rnacentral_pipeline.databases.data import Entry, Exon, SequenceRegion
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs


def accession(info):
    return "EXPRESSIONATLAS:" + info["GeneID"]


def primary_id(info):
    return "EXPRESSIONATLAS:" + info["GeneID"]


def taxid(info):
    taxid = info["taxid"][0]
    return int(taxid)


def species(info):
    return phy.species(info["taxid"][0])


def lineage(info):
    return phy.lineage(info["taxid"][0])


def common_name(info):
    return phy.common_name(info["taxid"][0])


def url(experiment):
    return "https://www.ebi.ac.uk/gxa/experiments/" + experiment


def region_builder(info):
    print(info["region_start"], info["region_stop"], info["strand"], info["urs_taxid"])
    return [
        SequenceRegion(
            chromosome=info["chromosome"][0],
            strand=info["strand"][0],
            exons=[
                Exon(start=start, stop=stop)
                for start, stop in zip(info["region_start"], info["region_stop"])
            ],
            assembly_id=info["assembly_id"][0],
            coordinate_system="1-start, fully-closed",
        )
    ]


def references(interactions):
    refs = set()
    for interaction in interactions:
        refs.update(interaction.publications)
    refs.add(pubs.reference(24234451))
    return list(refs)


def as_entry(info, experiment):
    print(info["seq"][0])
    print([""] if info["Gene Name"] == [None] else info["Gene Name"])
    return Entry(
        primary_id=primary_id(info),
        accession=accession(info),
        ncbi_tax_id=taxid(info),
        database="EXPRESSIONATLAS",
        sequence=info["seq"][0],
        regions=region_builder(info),
        rna_type=info["rna_type"][0],
        url=url(experiment),
        seq_version="1",
        description=info["description"][0],
        species=species(info),
        common_name=common_name(info),
        lineage=lineage(info),
        gene=info["GeneID"][0],
        gene_synonyms=[""]
        if info.get("Gene Name", [None]) == [None]
        else info["Gene Name"],
    )
