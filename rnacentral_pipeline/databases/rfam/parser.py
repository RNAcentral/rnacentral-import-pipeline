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

import csv
import logging
import tempfile
import typing as ty
from pathlib import Path

from Bio import SeqIO
from Bio.SeqIO import SeqRecord

from ..data import Entry
from . import helpers

LOGGER = logging.getLogger(__name__)


def as_entry(
    families: ty.Dict[str, ty.Dict[str, str]], sequences, data: ty.Dict[str, str]
) -> ty.Optional[Entry]:
    """
    Turn an entry from the JSON file into a Entry object for writing.
    """

    family = families[data["rfam_acc"]]
    try:
        sequence = helpers.sequence(sequences, data)
    except Exception:
        LOGGER.warn("Could not load sequence for %s", helpers.sequence_id(data))
        return None

    return Entry(
        primary_id=helpers.primary_id(data),
        accession=helpers.accession(data),
        ncbi_tax_id=helpers.taxid(data),
        database="RFAM",
        sequence=sequence,
        regions=[],
        rna_type=helpers.rna_type(family),
        url=helpers.url(family),
        note_data=helpers.note(data),
        species=helpers.species(data),
        lineage=helpers.lineage(data),
        optional_id=helpers.optional_id(family),
        product=helpers.product(family),
        parent_accession=helpers.parent_accession(data),
        project="RFAM",
        description=helpers.description(family, data),
        mol_type=helpers.mol_type(data),
        seq_version=helpers.seq_version(data),
        is_composite="N",
        location_start=helpers.location_start(data),
        location_end=helpers.location_end(data),
        references=helpers.references(data),
    )


def load_mapping(handle: ty.TextIO) -> ty.Dict[str, ty.Dict[str, str]]:
    reader = csv.DictReader(handle, delimiter="\t")
    data = {}
    for row in reader:
        data[row["id"]] = row
    return data


def dedup_sequences(path: Path) -> ty.Iterator[SeqRecord]:
    seen = set()
    for record in SeqIO.parse(str(path), "fasta"):
        if record.id in seen:
            LOGGER.warn("Duplicate id seen %s", record.id)
            continue
        seen.add(record.id)
        yield record


def parse(
    family_file: ty.TextIO, sequence_info: ty.TextIO, fasta: Path
) -> ty.Iterable[Entry]:
    """
    Parse the JSON file of Rfam data and produce a generator of all Entry
    objects in the file.
    """

    with tempfile.NamedTemporaryFile() as temp:
        SeqIO.write(dedup_sequences(fasta), temp.name, "fasta")
        temp.flush()
        indexed = SeqIO.index(temp.name, "fasta")
        families = load_mapping(family_file)
        reader = csv.DictReader(sequence_info, delimiter="\t")
        total = 0
        missing = 0
        for row in reader:
            total += 1
            entry = as_entry(families, indexed, row)
            if entry is None:
                missing += 1
                continue
            yield entry

    if missing:
        LOGGER.warn(
            "Did not load %i of %i sequences, yeilding: %i",
            missing,
            total,
            total - missing,
        )
    if float(missing) / float(total) >= 0.5:
        raise ValueError("Failed to find too many sequences")
