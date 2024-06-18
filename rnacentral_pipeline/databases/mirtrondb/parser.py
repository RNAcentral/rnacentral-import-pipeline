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
import re
import typing as ty

from rnacentral_pipeline.databases.data import Entry, Region, RelatedSequence

RNA_TYPES = {
    "mature": "SO:0000276",
    "precursor": "SO:0001244",
}

SPECIES = {
    "A. thaliana": 3702,
    "B. taurus": 9913,
    "C. elegans": 6239,
    "C. familiaris": -1,
    "D. melanogaster": 7227,
    "D. pseudoobscura": -1,
    "D. rerio": -1,
    "D. simulans": -1,
    "G. gallus": -1,
    "H. sapiens": 9606,
    "M. esculenta": -1,
    "M. mulatta": -1,
    "M. musculus": -1,
    "M. truncatula": -1,
    "O. sativa": -1,
    "P. troglodytes": -1,
    "S. Italica": 4555,
    "S. scrofa": 9823,
}


def text_value(row: ty.Dict[str, str], name: str) -> str | None:
    value = row[name].strip()
    if value == "-":
        return None
    return value


def regions(entry: ty.Dict[str, str]) -> ty.List[Region]:
    return None


def tax_id(entry: ty.Dict[str, str]) -> int:
    return SPECIES[entry["specie"].strip()]


def parse(handle: ty.IO):
    blank = handle.readline()
    assert not blank
    notification = handle.readline()
    assert notification == "##mirtronDB tabular format"
    reader = csv.DictReader(handle, delimiter="\t")

    pre = {}
    mature = {}
    for raw in reader:
        entry = Entry(
            primary_id=f"MIRTRONDB:{raw['id']}",
            accesion=raw["name"],
            ncbi_tax_id=tax_id(raw),
            database="MIRTRONDB",
            sequence=raw["sequence"],
            regions=regions(raw),
            rna_type=RNA_TYPES[raw["type"]],
            url=f"http://mirtrondb.cp.utfpr.edu.br/fetch_details.php?mrt_details={raw['name']}",
            seq_version="1",
            gene=raw["host gene"],
        )
        assert entry.accession not in pre
        assert entry.accession not in mature
        if raw["type"] == "mature":
            mature[entry.accession] = entry
        elif raw["type"] == "precursor":
            pre[entry.accession] = entry
        else:
            raise ValueError(f"Cannot handle {raw}")

    for id, entry in mature.items():
        pre_id = re.sub(r"-[35]p", "", id)
        if pre_id not in pre:
            continue
        pre_entry = pre[pre_id]
        pre_entry.related_sequences.append(
            RelatedSequence(
                sequence_id=id,
                relationship="mature_product",
                coordinates=find_coord(pre_entry.sequence, entry.sequence),
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
