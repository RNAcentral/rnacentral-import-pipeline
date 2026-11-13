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

"""
Map free text (e.g. a disease name) to an ontology term using the OLS4 search
API. Matching is deliberately conservative: manual override or exact label.
Exact synonym matching (the best ranked term when several share it) exists but
is off by default; see docs/disease-ontology-plan.md. Anything else is
unmapped.
"""

import asyncio
import csv
import logging
from collections import Counter
from typing import Dict, Optional, Tuple

from furl import furl

from . import fetch

LOGGER = logging.getLogger(__name__)

SEARCH_URL = "https://www.ebi.ac.uk/ols4/api/search"
HEADER = ["text", "ontology_term_id", "label", "match_type"]
# Diseases first; HPO catches conditions MONDO only has as phenotypes (sepsis,
# obesity, hypertension).
DEFAULT_ONTOLOGIES = "mondo,hp"
# The exact match is often not in the first 10 (e.g. "colon cancer").
ROWS = 50


def normalise(text: str) -> str:
    return " ".join(text.split()).casefold()


def load_overrides(path) -> Dict[str, str]:
    """
    Read the overrides TSV (header: text, ontology_term_id) into a dict keyed
    on the normalised text.
    """
    overrides = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            term_id = (row["ontology_term_id"] or "").strip()
            if not term_id:
                raise ValueError(f"Override for {row['text']!r} has no term id")
            overrides[normalise(row["text"])] = term_id
    return overrides


def _search(
    key: str, ontology: str, synonyms: bool = False
) -> Optional[Tuple[str, str, str]]:
    """
    One OLS search in one ontology. Docs come back best match first, and that
    order is the tie-break when several terms share an exact synonym.
    """
    url = furl(SEARCH_URL).add(
        {
            "q": key,
            "ontology": ontology.lower(),
            "rows": ROWS,
            "queryFields": "label,synonym",
            "fieldList": "obo_id,label,synonym,is_defining_ontology",
        }
    )
    result = asyncio.run(fetch.query_ols(url))
    prefix = ontology.upper() + ":"
    docs = [
        doc
        for doc in result["response"]["docs"]
        if (doc.get("obo_id") or "").startswith(prefix)
    ]

    for doc in docs:
        if normalise(doc.get("label") or "") == key:
            return (doc["obo_id"], doc["label"], "label_exact")
    for doc in docs if synonyms else []:
        if any(normalise(syn) == key for syn in doc.get("synonym") or []):
            return (doc["obo_id"], doc.get("label"), "synonym_exact")
    return None


def map_text(
    text: str,
    ontology: str = DEFAULT_ONTOLOGIES,
    overrides: Optional[Dict[str, str]] = None,
    synonyms: bool = False,
) -> Tuple[Optional[str], Optional[str], str]:
    """
    Returns (ontology_term_id, label, match_type). `ontology` is a comma
    separated list tried in order, so later ones are only a fallback.
    `overrides` is keyed on the normalised text. An override has no label,
    since OLS is not queried.
    """
    key = normalise(text)
    if overrides and key in overrides:
        return (overrides[key], None, "override")
    for name in ontology.split(","):
        found = _search(key, name.strip(), synonyms)
        if found:
            return found
    return (None, None, "unmapped")


def map_file(
    handle, output, ontology=DEFAULT_ONTOLOGIES, overrides_path=None, synonyms=False
):
    """
    Read one string per line from `handle`, dedupe on the normalised form and
    write the mapping TSV to the `output` handle.
    """
    overrides = load_overrides(overrides_path) if overrides_path else None
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADER)
    seen = set()
    counts = Counter()
    for line in handle:
        text = line.rstrip("\r\n")
        key = normalise(text)
        if not key or key in seen:
            continue
        seen.add(key)
        term_id, label, match_type = map_text(text, ontology, overrides, synonyms)
        counts[match_type] += 1
        writer.writerow([text, term_id or "", label or "", match_type])
    LOGGER.info("Mapped %i strings: %s", len(seen), dict(counts))
