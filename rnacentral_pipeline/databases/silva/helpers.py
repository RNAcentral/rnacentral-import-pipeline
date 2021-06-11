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
import logging
import typing as ty
import operator as op

from rnacentral_pipeline.databases import data
import rnacentral_pipeline.databases.helpers.publications as pubs
import rnacentral_pipeline.databases.helpers.phylogeny as phy

url = op.itemgetter("silvaUri")

CLASS_PATTERN = re.compile(r";\s*$")

KNOWN_TYPES = {
    "rRNA": "SO:0000252",
    "rRNA_12S": "SO:0002128",
    "rRNA_16S": "SO:0001000",
    "rRNA_18S": "SO:0000407",
    "small_subunit_rRNA": "SO:0000650",
    "large_subunit_rRNA": "SO:0000651",
    "rRNA_23S": "SO:0001001",
    "rRNA_26S": "SO:0000653",  # TODO: Get a better term
    "rRNA_28S": "SO:0000653",
}

RRNA_NAME_MAPPING = {
    "rRNA": "rRNA",
    "rRNA_12S": "mitochondrial SSU rRNA",
    "rRNA_16S": "bacterial SSU rRNA",
    "rRNA_18S": "SSU rRNA",
    "small_subunit_rRNA": "SSU rRNA",
    "large_subunit_rRNA": "LSU rRNA",
    "rRNA_23S": "bacterial LSU rRNA",
    "rRNA_26S": "eukaryotic LSU rRNA",
    "rRNA_28S": "eukaryotic LSU rRNA",
}

LOGGER = logging.getLogger(__name__)


def primary_id(row) -> str:
    return "SILVA:%s:%s" % (row["insdcAccession"], row["location"])


def taxid(row) -> int:
    return int(row["ncbiTaxId"])


def sequence(row) -> str:
    return row["sequence"].replace("U", "T")


def inference(row):
    value = row["classification"]
    value = value.replace(";", "; ")
    return re.sub(CLASS_PATTERN, "", value)


def version(row) -> str:
    _, version = row["insdcAccession"].split(".", 1)
    return version


def rna_type(row) -> str:
    given = row["type"]
    if given in KNOWN_TYPES:
        return KNOWN_TYPES[given]
    raise ValueError("Unknown RNA type")


def lineage(taxonomy, row) -> str:
    tid = taxid(row)
    if tid in taxonomy:
        return taxonomy[tid].lineage
    return phy.lineage(tid)


def species(taxonomy, row) -> str:
    tid = taxid(row)
    if tid in taxonomy:
        return taxonomy[tid].name
    return phy.species(tid)


def description(taxonomy, row) -> str:
    organism = species(taxonomy, row)
    rrna = RRNA_NAME_MAPPING[row["type"]]
    return f"{organism} {rrna}"


def as_entry(taxonomy, row) -> ty.Optional[data.Entry]:
    try:
        return data.Entry(
            primary_id=primary_id(row),
            accession=primary_id(row),
            ncbi_tax_id=taxid(row),
            database="SILVA",
            sequence=sequence(row),
            regions=[],
            rna_type=rna_type(row),
            url=url(row),
            seq_version=version(row),
            species=species(taxonomy, row),
            lineage=lineage(taxonomy, row),
            references=[
                pubs.reference("doi:10.1093/nar/gks1219"),
            ],
            inference=inference(row),
            description=description(taxonomy, row),
        )
    except phy.FailedTaxonId:
        LOGGER.warning("Could not get phylogeny info for %s", row)
        return None
    except phy.UnknownTaxonId:
        return None
