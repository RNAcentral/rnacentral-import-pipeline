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

import csv
import operator as op
import re
import typing as ty

from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source

DOMAINS = {
    "arch": ("A", 2157),
    "euk": ("E", 2759),
    "bact": ("B", 2),
}

TYPES = {
    "Ala": "SO:0000254",
    "Arg": "SO:0001036",
    "Asn": "SO:0000256",
    "Asp": "SO:0000257",
    "Cys": "SO:0000258",
    "Gln": "SO:0000259",
    "Glu": "SO:0000260",
    "Gly": "SO:0000261",
    "His": "SO:0000262",
    "Ile": "SO:0000263",
    "Ile2": "SO:0000263",
    "Leu": "SO:0000264",
    "Lys": "SO:0000265",
    "Met": "SO:0000266",
    "Phe": "SO:0000267",
    "Pro": "SO:0000268",
    "SeC": "SO:0005857",
    "Ser": "SO:0000269",
    "Thr": "SO:0000270",
    "Trp": "SO:0000271",
    "Tyr": "SO:0000272",
    "Val": "SO:0000273",
    "fMet": "SO:0000266",  # TODO: Improve
    "iMet": "SO:0000266",  # TODO: Improve
}


def parse_model(handle) -> ModelInfo:
    name: ty.Optional[str] = None
    length: ty.Optional[str] = None
    for line in handle:
        line = line.strip()
        if line == "CM":
            break
        key, value = re.split("\s+", line, maxsplit=1)
        if key == "NAME":
            name = value
        if key == "CLEN":
            length = value

    # TODO: Handle parsing organelle models
    loc = "cellular"

    if not name:
        raise ValueError("Invalid name")

    if not length:
        raise ValueError("Invalid length for: %s" % name)

    domain, trna_type = name.split("-")
    if domain not in DOMAINS:
        raise ValueError("Cannot find taxid for model: " + name)
    if trna_type not in TYPES:
        raise ValueError("Cannot find SO term for model: " + name)

    short_domain, taxid = DOMAINS[domain]
    so_term = TYPES[trna_type]
    model_name = "%s-%s" % (short_domain, trna_type)

    return ModelInfo(
        model_name=model_name,
        so_rna_type=so_term,
        taxid=taxid,
        source=Source.gtrnadb,
        length=int(length),
        cell_location=loc,
        basepairs=None,
    )


def parse(handle, extra=None):
    for line in handle:
        if line.startswith("INFERNAL"):
            yield parse_model(handle)


def write(handle, output):

    data = parse(handle)
    data = map(op.methodcaller("writeable"), data)
    csv.writer(output).writerows(data)
