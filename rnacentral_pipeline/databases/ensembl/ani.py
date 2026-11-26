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

import json
import time
import typing as ty
from urllib.request import urlopen

API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/dataset_report?page_size={}"
BATCH = 100


def assemblies(data: dict) -> ty.List[str]:
    """Every assembly in species.json, so a species can fall back past a failed one."""
    return sorted(
        {acc for info in data["species"].values() for acc in info.get("assemblies", {})}
    )


def reports(accessions: ty.List[str]) -> ty.Iterable[dict]:
    for start in range(0, len(accessions), BATCH):
        batch = accessions[start : start + BATCH]
        with urlopen(API.format(",".join(batch), BATCH), timeout=120) as handle:
            yield from json.load(handle).get("reports", [])
        # NCBI allows 5 requests a second without an API key
        time.sleep(0.4)


SPECIES_ANI = 95.0


def failed(found: ty.Iterable[dict]) -> ty.Iterable[str]:
    """Genomes that are not the species named: Failed, or closer to another species."""
    for report in found:
        ani = report.get("average_nucleotide_identity") or {}
        own = (ani.get("submitted_ani_match") or {}).get("ani")
        best = (ani.get("best_ani_match") or {}).get("ani")
        # NCBI calls Lacrimispora indolis Inconclusive at 82.76% own, 99.99% S. nepalensis
        mixed_up = own is not None and own < SPECIES_ANI <= (best or 0)
        if ani.get("taxonomy_check_status") == "Failed" or mixed_up:
            yield report["accession"]
