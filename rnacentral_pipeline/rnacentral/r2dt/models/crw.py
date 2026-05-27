# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source

SO_TERM_MAPPING = {
    "rRNA_16S": "SO:0000650",
    "rRNA_5S": "SO:0000652",
    "group_I_intron": "SO:0000587",
    "group_II_intron": "SO:0000603",
    "large_subunit_rRNA": "SO:0000651",
    "small_subunit_rRNA": "SO:0000650",
    "mt_LSU_rRNA": "SO:0002345",
    "mt_SSU_rRNA": "SO:0002344",
    "rRNA_18S": "SO:0000407",
    "rRNA_21S": "SO:0002345",
    "rRNA_23S": "SO:0001001",
}


def load_stats(handle) -> dict:
    stats = {}
    if handle is None:
        return stats
    for row in csv.reader(handle):
        if row:
            stats[row[0]] = {"length": int(row[1]), "basepairs": int(row[2])}
    return stats


def as_so_term(raw):
    if raw in SO_TERM_MAPPING:
        return SO_TERM_MAPPING[raw]
    raise ValueError("Unknown RNA type: " + raw)


def as_taxid(raw):
    if raw == "501083":
        return 126
    if raw in {"600001", "600002", "600003"}:
        return 562
    if raw in {"600101", "600102"}:
        return 2238
    if raw in {"600301", "600302"}:
        return 4932
    if raw in {"600201", "600202"}:
        return 274
    return int(raw)


def parse(handle, extra=None):
    stats = load_stats(extra)
    for row in csv.DictReader(handle, delimiter="\t"):
        model_name = row["model_name"]
        so_type_name = row["rna_type"]
        if so_type_name == "mt_rRNA":
            so_type_name = "mt_SSU_rRNA" if ".16." in model_name else "mt_LSU_rRNA"
        try:
            so_term_val = as_so_term(so_type_name)
        except ValueError:
            continue
        model_stats = stats.get(model_name, {})
        yield ModelInfo(
            model_name=model_name,
            so_rna_type=so_term_val,
            taxid=as_taxid(row["taxid"]),
            source=Source.crw,
            length=model_stats.get("length"),
            basepairs=model_stats.get("basepairs"),
            cell_location=None,
        )
