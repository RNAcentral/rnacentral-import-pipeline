# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import csv
import logging
import re
import typing as ty

from rnacentral_pipeline.databases.data import Entry, RelatedCoordinate, RelatedSequence

LOGGER = logging.getLogger(__name__)

RNA_TYPES = {
    "mature": "SO:0000276",
    "precursor": "SO:0001244",
}

SPECIES = {
    "A. thaliana": 3702,
    "B. taurus": 9913,
    "C. elegans": 6239,
    "C. familiaris": 9615,
    "D. melanogaster": 7227,
    "D. pseudoobscura": 7237,
    "D. rerio": 7955,
    "D. simulans": 7240,
    "G. gallus": 9031,
    "H. sapiens": 9606,
    "M. esculenta": 3983,
    "M. mulatta": 9544,
    "M. musculus": 10090,
    "M. truncatula": 3880,
    "O. sativa": 4530,
    "P. troglodytes": 9598,
    "S. Italica": 4555,
    "S. scrofa": 9823,
}


def text_value(row: ty.Dict[str, str], name: str) -> str | None:
    value = row[name].strip()
    if value == "-":
        return None
    return value


def find_coords(id: str, target: str, query: str) -> ty.List[RelatedCoordinate]:
    if query not in target:
        LOGGER.warn(f"Mature not found in precusor for %s", id)
        return []
    start = target.index(query)
    return [RelatedCoordinate(start=start, stop=start + len(query))]


def parse(handle: ty.IO):
    blank = handle.readline().strip()
    assert not blank, f"Invalid first line `{blank}`"
    notification = handle.readline().strip()
    assert notification == "##mirtronDB tabular format"
    reader = csv.DictReader(handle, delimiter="\t")

    pre = {}
    mature = {}
    for raw in reader:
        rna_type = raw["type"].strip()
        name = raw["name"].strip()
        species = raw["specie"].strip()
        description = f"{species} {name} {rna_type} miRNA"
        if raw["host gene"].strip():
            description += f" ({raw['host gene'].strip()})"
        entry = Entry(
            primary_id=f"MIRTRONDB:{raw['id'].strip()}",
            accession=name,
            ncbi_tax_id=SPECIES[species],
            database="MIRTRONDB",
            sequence=raw["sequence"].strip(),
            regions=[],
            rna_type=RNA_TYPES[rna_type],
            url=f"http://mirtrondb.cp.utfpr.edu.br/fetch_details.php?mrt_details={name}",
            seq_version="1",
            gene=raw["host gene"].strip(),
            description=description,
        )
        assert entry.accession not in pre
        assert entry.accession not in mature

        if rna_type == "mature":
            mature[entry.accession] = entry
        elif rna_type == "precursor":
            pre[entry.accession] = entry
        else:
            raise ValueError(f"Cannot handle {raw}")

    for id, entry in mature.items():
        pre_id = re.sub(r"-[35]p", "", id)
        if pre_id not in pre:
            LOGGER.warn("Failed to find precursor for %s", id)
            continue
        pre_entry = pre[pre_id]
        pre_entry.related_sequences.append(
            RelatedSequence(
                sequence_id=id,
                relationship="mature_product",
                coordinates=find_coords(id, pre_entry.sequence, entry.sequence),
            )
        )

        entry.related_sequences.append(
            RelatedSequence(
                sequence_id=pre_id,
                relationship="precursor",
            )
        )

    yield from pre.values()
    yield from mature.values()
