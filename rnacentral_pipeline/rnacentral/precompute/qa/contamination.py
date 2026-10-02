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

import re
import typing as ty

from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult

ALLOWED_FAMILIES = {
    "RF00177",  # Bacterial small subunit ribosomal RNA
    "RF02541",  # Bacterial large subunit ribosomal RNA
    "RF01959",  # Archaeal small subunit ribosomal RNA
    "RF02540",  # Archaeal large subunit ribosomal RNA
}

GENERIC_DOMAINS = {
    "unclassified sequences",
    "artificial sequences",
    "miscellaneous sequences",
    "other sequences",
}


def is_ignorable_organelle_conflict(rna_type: str, data: Sequence) -> bool:
    """
    This can ignore any conflict where the sequence probably comes from a
    mitochondrion or chloroplast but matches a bacterial rRNA. In that case
    we do not warn since this is expected from evolution (both organelles
    have a bacterial evolutionary origin).
    """
    if rna_type != "rRNA" or data.rfam_hits[0].model not in ALLOWED_FAMILIES:
        return False
    return data.is_mitochondrial() or data.is_chloroplast()


def is_generic_domain(domains: ty.Set[str]) -> bool:
    """
    Check if any domain in the given set is a generic domain (ie
    unclassified sequences, etc).
    """
    return bool(domains & GENERIC_DOMAINS)


def message(data: Sequence) -> str:
    """
    Produce a warning message about the issue detected. This assumes that
    there was a warning.
    """

    common_name = {acc.common_name for acc in data.accessions}
    common_name = {c.lower() for c in common_name if c}

    # No common_name or species to identify the sequence with (e.g. no active
    # accessions carry either) - leave it blank, the final whitespace cleanup
    # below collapses "This  sequence" down to "This sequence".
    sequence_name = ""
    if len(common_name) == 1:
        sequence_name = common_name.pop()
    else:
        species = {acc.species for acc in data.accessions if acc.species}
        if species:
            sequence_name = ", ".join(sorted(species))
            sequence_name = f"<i>{sequence_name}</i>"

    model_domain = data.rfam_hits[0].model_domain
    model_url = data.rfam_hits[0].url
    model_name = data.rfam_hits[0].model_name

    msg = (
        "This {sequence_name} sequence matches a {match_domain} "
        'Rfam model (<a href="{model_url}">{model_name}</a>). '
        '<a href="{help_url}">Learn more &rarr;</a>'
    ).format(
        sequence_name=sequence_name,
        match_domain=model_domain,
        model_url=model_url,
        model_name=model_name,
        help_url="/help/qc",
    )

    return re.sub(r"\s+", " ", msg)


def validate(rna_type: str, sequence: Sequence) -> QaResult:
    if not sequence.has_unique_hit():
        return QaResult.ok("possible_contamination")

    hit = sequence.rfam_hits[0]
    domains = sequence.domains()
    if not hit.model_domain or not domains or is_generic_domain(domains):
        return QaResult.ok("possible_contamination")

    if hit.model_domain not in domains and not is_ignorable_organelle_conflict(
        rna_type, sequence
    ):
        return QaResult.not_ok("possible_contamination", message(sequence))
    return QaResult.ok("possible_contamination")
