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

import collections as coll
import logging
import re
import typing as ty

from rnacentral_pipeline.databases.data import IdReference
from rnacentral_pipeline.databases.helpers.publications import reference

LOGGER = logging.getLogger(__name__)


def primary_id(data: ty.Dict[str, str]) -> str:
    return data["rfam_acc"]


def accession(data: ty.Dict[str, str]) -> str:
    return f"{data['rfamseq_acc']}.{data['version']}:{data['seq_start']}..{data['seq_end']}:rfam"


def taxid(data: ty.Dict[str, str]) -> int:
    return int(data["ncbi_id"])


def sequence_id(data: ty.Dict[str, str]) -> str:
    return f"{data['rfamseq_acc']}/{data['seq_start']}-{data['seq_end']}"


def sequence(sequences, data: ty.Dict[str, str]) -> str:
    seq_id = sequence_id(data)
    # print((seq_id, data["sequence_type"]))
    sequence = str(sequences[seq_id].seq)
    return sequence.upper().replace("U", "T")


def seq_version(data: ty.Dict[str, str]) -> str:
    return data["version"]


def rna_type(family: ty.Dict[str, str]) -> str:
    so_terms = family["so_terms"]
    if "," in so_terms:
        so_terms = so_terms.split(",")[0]
    if "," in so_terms:
        so_terms = so_terms.split(",")[0]
    assert re.match(r"^SO:\d+$", so_terms)
    return so_terms


def url(family: ty.Dict[str, str]) -> str:
    return f"http://rfam.org/family/{family['id']}"


def species(data: ty.Dict[str, str]) -> str:
    return data["species"]


def lineage(data):
    base = re.sub(r"\.$", "", data["tax_string"])
    return "; ".join(base.split(" ") + [data["species"]])


def references(data: ty.Dict[str, str]) -> ty.List[IdReference]:
    refs = [reference(29112718)]
    pmids = data.get("PMIDS", "").split(",")
    for pmid in pmids:
        try:
            refs.append(reference(pmid))
        except:
            LOGGER.warn("Failed to parse reference '%s'", pmid)
    return refs


def note(data: ty.Dict[str, str]):
    result = coll.defaultdict(list)
    result["Alignment"] = data["sequence_type"]
    for xref in data["dbxrefs"].split(","):
        db, _ = xref.split(":", 1)
        result[db].append(xref)
    return result


def experiment(data) -> str:
    return " ".join(p.external_id for p in references(data))


def description(family, data: ty.Dict[str, str]) -> str:
    return f"{species(data)} {product(family)}"


def optional_id(family: ty.Dict[str, str]) -> str:
    return family["name"]


def product(family: ty.Dict[str, str]) -> str:
    return family["pretty_name"]


def parent_accession(data: ty.Dict[str, str]) -> str:
    return data["rfamseq_acc"]


def mol_type(data: ty.Dict[str, str]) -> str:
    return data["sequence_type"]


def location_start(data: ty.Dict[str, str]) -> int:
    return int(data["seq_start"])


def location_end(data: ty.Dict[str, str]) -> int:
    return int(data["seq_end"])
