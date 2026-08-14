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

import logging
import typing as ty
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pub
from rnacentral_pipeline.databases.helpers import r2dt
from rnacentral_pipeline.rnacentral.r2dt.data import SO_RNA_NAME_LOOKUP

LOGGER = logging.getLogger(__name__)


def primary_id(info) -> str:
    return "CRW:" + info.model_name


def sequence(info, sequences) -> str:
    return str(sequences[info.model_name].seq).upper().replace("U", "T")


def description(info) -> str:
    # The lookup holds internal names ("cytosolic_SSU_rRNA"); descriptions
    # are user facing, so space them out.
    rna_type = SO_RNA_NAME_LOOKUP[info.so_rna_type].replace("_", " ")
    return f"{phy.species(info.taxid)} {rna_type}"


def as_entry(info, sequences) -> ty.Optional[data.Entry]:
    """
    Build an entry from the ModelInfo r2dt.models.crw parses out of R2DT's
    crw-metadata.tsv. No organelle: the metadata has no such column and
    r2dt_models never stored one either.
    """
    try:
        return data.Entry(
            primary_id=primary_id(info),
            accession=primary_id(info),
            ncbi_tax_id=info.taxid,
            database="CRW",
            regions=[],
            rna_type=info.so_rna_type,
            sequence=sequence(info, sequences),
            url="",
            seq_version="1",
            description=description(info),
            species=phy.species(info.taxid),
            common_name=phy.common_name(info.taxid),
            lineage=phy.lineage(info.taxid),
            references=[
                pub.reference(11869452),
            ],
        )
    except Exception as err:
        LOGGER.warning("Could not generate entry for %s", info.model_name)
        LOGGER.exception(err)
        return None


def fasta_entries(directory: Path) -> ty.Iterable[SeqRecord]:
    return r2dt.fasta_entries([directory])
