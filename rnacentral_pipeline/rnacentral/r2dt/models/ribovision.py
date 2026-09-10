# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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
import typing as ty

from rnacentral_pipeline.databases.helpers.phylogeny import FailedTaxonId, taxid
from rnacentral_pipeline.rnacentral.r2dt.data import ModelInfo, Source


def load_stats(handle) -> dict:
    stats = {}
    if handle is None:
        return stats
    for row in csv.reader(handle):
        if row:
            stats[row[0]] = {"length": int(row[1]), "basepairs": int(row[2])}
    return stats


def lookup_taxid(species):
    if species == "Deinococcus radiodurans":
        return 1299
    if species == "Bacillus subtilis":
        return 1423
    if species == "Spinacia oleracea":
        return 3562
    if species == "Dictyostelium discoideum":
        return 44689
    if species == "Drosophilla melanogaster":
        return 7227
    if species == "Escherichia coli":
        return 562
    if species == "Haloarcula marismortui":
        return 2238
    if species == "Homo sapiens":
        return 9606
    if species == "Kluyveromyces lactis":
        return 28985
    if species == "Mycolicibacterium smegmatis":
        return 1772
    if species == "Mycobacterium tuberculosis":
        return 1773
    if species == "Plasmodium falciparum":
        return 5833
    if species == "Staphylococcus aureus":
        return 1280
    if species == "Saccharomyces cerevisiae":
        return 4932
    if species == "Tetrahymena thermophila":
        return 5911
    if species == "Thermus thermophilus":
        return 274
    if species == "Drosophila melanogaster":
        return 7227
    if species == "Trypanosoma brucei":
        return 5691
    try:
        return taxid(species)
    except FailedTaxonId:
        raise ValueError("Unknown species name: " + species)


def as_location(raw):
    if not raw:
        return None
    if raw == "mito":
        return "Mitochondrion"
    if raw == "chloroplast":
        return "Chloroplast"
    raise ValueError("Unknown raw location: " + raw)


def location_of(row):
    """
    The location sits in a fourth column that metadata.tsv has no header for,
    so DictReader collects it under None. Rows without one omit it entirely.
    """
    extra = row.get(None) or []
    return as_location(extra[0] if extra else None)


def so_term(row):
    if row["model_name"] == "DR_LSU_3D":
        return "SO:0001001"
    if (
        "LSU" in row["model_name"]
        or "23S" in row["model_name"]
        or "28S" in row["model_name"]
        or row["model_name"] == "F3H4_G_18380"
        or row["model_name"] == "GC14_75"
    ):
        return "SO:0000651"
    if "SSU" in row["model_name"]:
        return "SO:0000650"
    raise ValueError("Could not figure out SO term for: %s" % row)


def parse(handle, extra=None):
    stats = load_stats(extra)
    for row in csv.DictReader(handle, delimiter="\t"):
        so_term_id = so_term(row)

        if not row.get("taxid", None):
            tid = lookup_taxid(row["species"])
        else:
            tid = int(row["taxid"])

        model_stats = stats.get(row["model_name"], {})
        yield ModelInfo(
            model_name=row["model_name"],
            so_rna_type=so_term_id,
            taxid=tid,
            source=Source.ribovision,
            cell_location=location_of(row),
            length=model_stats.get("length"),
            basepairs=model_stats.get("basepairs"),
        )
