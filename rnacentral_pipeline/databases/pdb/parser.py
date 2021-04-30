# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import logging
import typing as ty

from rnacentral_pipeline.databases import data

from rnacentral_pipeline.databases.pdb.data import ChainInfo
from rnacentral_pipeline.databases.pdb.data import ReferenceMapping
from rnacentral_pipeline.databases.pdb import helpers

LOGGER = logging.getLogger(__name__)


def as_entry(info: ChainInfo, reference_mapping: ReferenceMapping):
    return data.Entry(
        primary_id=info.pdb_id.upper(),
        accession=info.accession(),
        ncbi_tax_id=helpers.taxid(info),
        database="PDBE",
        sequence=helpers.sequence(info),
        regions=[],
        rna_type=helpers.rna_type(info),
        url=helpers.url(info),
        seq_version="1",
        note_data=helpers.note_data(info),
        product=helpers.product(info),
        optional_id=info.chain_id,
        description=helpers.description(info),
        species=helpers.species(info),
        lineage=helpers.lineage(info),
        parent_accession=info.pdb_id.upper(),
        references=helpers.references_for(info, reference_mapping),
    )


def parse(
    rna_chains: ty.List[ChainInfo],
    reference_mapping: ReferenceMapping,
) -> ty.Iterator[data.Entry]:
    disqualified = {"mRNA": 0, "other": 0}
    for chain in rna_chains:
        if helpers.is_mrna(chain):
            LOGGER.debug("Disqualifing %s", chain)
            disqualified["mRNA"] += 1
            continue

        if not helpers.is_ncrna(chain):
            LOGGER.debug("Skipping %s", chain)
            disqualified["other"] += 1
            continue

        try:
            yield as_entry(chain, reference_mapping)
        except helpers.InvalidSequence:
            LOGGER.warn(f"Invalid sequence for {chain}")
    LOGGER.info("Disqualified %i mRNA chains", disqualified["mRNA"])
    LOGGER.info("Disqualified %i non ncRNA chains", disqualified["other"])
