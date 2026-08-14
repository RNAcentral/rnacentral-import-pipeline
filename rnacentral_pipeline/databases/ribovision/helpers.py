# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

LOGGER = logging.getLogger(__name__)

URL = "https://ribovision2.chemistry.gatech.edu/"

ORGANELLE_MAPPING = {
    "Mitochondrion": "mitochondria",
    "Cyanelle": "cyanelle",
    "Chloroplast": "chloroplast",
}


def primary_id(row: ty.Dict[str, ty.Any]) -> str:
    return "RIBOVISION:" + row["model_name"]


def taxid(row: ty.Dict[str, ty.Any]) -> int:
    return row["taxid"]


def sequence(row: ty.Dict[str, ty.Any], sequences) -> str:
    return str(sequences[row["model_name"]].seq).upper().replace("U", "T")


def organelle(row: ty.Dict[str, ty.Any]) -> ty.Optional[str]:
    cellular_location = row.get("cellular_location", None)
    if cellular_location is None:
        return None
    return ORGANELLE_MAPPING.get(cellular_location, None)


def description(row: ty.Dict[str, ty.Any]) -> str:
    name = phy.species(taxid(row))
    loc = organelle(row)
    rna_type = row["rna_type"]
    if loc:
        return f"{name} {loc} {rna_type}"
    return f"{name} {rna_type}"


def as_entry(row: ty.Dict[str, ty.Any], sequences) -> ty.Optional[data.Entry]:
    try:
        return data.Entry(
            primary_id=primary_id(row),
            accession=primary_id(row),
            ncbi_tax_id=taxid(row),
            database="RIBOVISION",
            regions=[],
            rna_type=row["so_term_id"],
            sequence=sequence(row, sequences),
            url=URL,
            seq_version="1",
            description=description(row),
            species=phy.species(taxid(row)),
            common_name=phy.common_name(taxid(row)),
            lineage=phy.lineage(taxid(row)),
            references=[
                pub.reference(39237196),
            ],
            organelle=organelle(row),
        )
    except Exception as err:
        LOGGER.warning("Could not generate entry for %s", row)
        LOGGER.exception(err)
        return None


def fasta_entries(directory: Path) -> ty.Iterable[SeqRecord]:
    """
    RiboVision models are split across an LSU and an SSU directory, each keeping
    its sequences alongside the bpseq they were generated from.
    """
    return r2dt.fasta_entries(sorted(directory.glob("ribovision-*/bpseq")))
