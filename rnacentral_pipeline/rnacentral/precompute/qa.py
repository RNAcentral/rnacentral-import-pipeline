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

INCOMPLETE_FAMILIES = {
    'RF00001',  # 5S ribosomal RNA
    'RF00002',  # 5.8S ribosomal RNA
    'RF00005',  # tRNA
    'RF00177',  # Bacterial small subunit ribosomal RNA
    'RF01959',  # Archaeal small subunit ribosomal RNA
    'RF01960',  # Eukaryotic small subunit ribosomal RNA
    'RF02540',  # Archaeal large subunit ribosomal RNA
    'RF02541',  # Bacterial large subunit ribosomal RNA
    'RF02542',  # Microsporidia small subunit ribosomal RNA
    'RF02543',  # Eukaryotic large subunit ribosomal RNA
}


EXPECTED_MATCHES = {
    'rRNA': {
        'RF00001',  # 5S ribosomal RNA
        'RF00002',  # 5.8S ribosomal RNA
        'RF00177',  # Bacterial small subunit ribosomal RNA
        'RF01959',  # Archaeal small subunit ribosomal RNA
        'RF01960',  # Eukaryotic small subunit ribosomal RNA
        'RF02540',  # Archaeal large subunit ribosomal RNA
        'RF02541',  # Bacterial large subunit ribosomal RNA
        'RF02542',  # Microsporidia small subunit ribosomal RNA
        'RF02543',  # Eukaryotic large subunit ribosomal RNA
        'RF02547',  # mito 5S RNA
    },
    'tRNA': {
        'RF00005',  # tRNA
        'RF01852',  # Selenocysteine tRNA
    },
}


def is_ignorable_mito_missing_trna(rna_type, data):
    """
    Check if the sequence is likely a mitochondrial tRNA so missing hits can be
    ignored in some cases.
    """

    if 'tRNA' == rna_type:
        return False

    return data.is_mitochondrial() and not data.hits


def is_ignorable_mito_conflict(rna_type, data):
    """
    This can ignore any conflict where the sequence probably comes from a
    mitochondria but it matches a bacterial rRNA. In that case we do not
    warn since this is expected from evolution.
    """

    return data.is_mitochondrial() and \
        'rRNA' == rna_type and \
        data.hits[0].rfam_model in {
            'RF00177',  # Bacterial small subunit ribosomal RNA
            'RF02541',  # Bacterial large subunit ribosomal RNA
            'RF01959',  # Archaeal small subunit ribosomal RNA
            'RF02540',  # Archaeal large subunit ribosomal RNA
        }


def possible_contamination(rna_type, data):
    """
    Check if the given sequence with the given RNA type is likely
    contamination. This checks if the domain of the sequence and the domain of
    the hits disagree.
    """

    # if not hits or len(hits) > 1:
    if not data.has_unique_hit():
        return False

    hit = data.hits[0]
    if not hit.model_domain:
        return False

    if not data.domains:
        return False

    return hit.model_domain not in data.domains and \
        not is_ignorable_mito_conflict(rna_type, data)


def incomplete_sequence(_, data):
    """
    Detect if the given sequence is incomplete. This will work by checking if
    the Rfam hits cover enough of the sequence and the hit is complete enough.
    """

    hits = data.hits
    if len(hits) != 1:
        return False

    if hits[0].rfam_model not in INCOMPLETE_FAMILIES:
        return False

    return hits[0].model_completeness <= 0.5 and \
        hits[0].sequence_completeness >= 0.9


def missing_rfam_match(rna_type, data):
    """
    Detect if the given sequence, with the given RNA type is missing an Rfam
    match.
    """

    rna_types = data.rna_types
    if len(rna_types) > 1:
        return False

    rna_type = rna_types.pop()
    if rna_type not in EXPECTED_MATCHES:
        return False

    if is_ignorable_mito_missing_trna(rna_type, data):
        return False

    required = EXPECTED_MATCHES[rna_type]
    hits = {h.rfam_model_id for h in data.hits}
    return not hits.intersection(required)
