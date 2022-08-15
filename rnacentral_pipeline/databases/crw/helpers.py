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
from pathlib import Path
import typing as ty

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pub

LOGGER = logging.getLogger(__name__)

ORGANELLE_MAPPING = {
    "Mitochondrion": "mitochondria",
    "Cyanelle": "cyanelle",
    "Chloroplast": "chloroplast",
}


def primary_id(row: ty.Dict[str, ty.Any]) -> str:
    return "CRW:" + row["model_name"]


def taxid(row: ty.Dict[str, ty.Any]) -> int:
    return row["taxid"]


def species(row: ty.Dict[str, ty.Any]) -> str:
    return phy.species(taxid(row))


def common_name(row: ty.Dict[str, ty.Any]) -> str:
    return phy.common_name(taxid(row))


def lineage(row: ty.Dict[str, ty.Any]) -> str:
    return phy.lineage(taxid(row))


def sequence(row: ty.Dict[str, ty.Any], sequences: ty.Dict[str, SeqRecord]) -> str:
    return str(sequences[row["model_name"]].seq).upper().replace("U", "T")


def description(row: ty.Dict[str, ty.Any]) -> str:
    name = species(row)
    loc = organelle(row)
    rna_type = row["rna_type"]
    if loc:
        return f"{name} {loc} {rna_type}"
    return f"{name} {rna_type}"


def organelle(row: ty.Dict[str, ty.Any]) -> ty.Optional[str]:
    cellular_location = row.get("cellular_location", None)
    if cellular_location is not None:
        return ORGANELLE_MAPPING.get(row["cellular_location"], None)
    else:
        return None

def as_entry(row: ty.Dict[str, ty.Any], sequences) -> ty.Optional[data.Entry]:
    try:
        return data.Entry(
            primary_id=primary_id(row),
            accession=primary_id(row),
            ncbi_tax_id=taxid(row),
            database="CRW",
            regions=[],
            rna_type=row["so_term_id"],
            sequence=sequence(row, sequences),
            url="",
            seq_version="1",
            description=description(row),
            species=species(row),
            common_name=common_name(row),
            lineage=lineage(row),
            references=[
                pub.reference(11869452),
            ],
            organelle=organelle(row),
        )
    except Exception as err:
        LOGGER.warn("Could not generate entry for %s", row)
        LOGGER.exception(err)
        return None


def fasta_entries(directory: Path) -> ty.Iterable[SeqRecord]:
    model_pattern = re.compile("crw-bpseq/(.+).bpseq ")
    for fasta_file in directory.glob("*.fasta"):
        with fasta_file.open("r") as raw:
            header, sequence, _ = raw.readlines()
            matches = re.search(model_pattern, header)
            if matches is None:
                raise ValueError(f"Could not get model id from {header}")
        yield SeqRecord(Seq(sequence.strip()), id=matches.group(1))
