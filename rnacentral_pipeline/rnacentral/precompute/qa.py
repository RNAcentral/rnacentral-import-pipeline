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


def is_ignorable_mito_conflict(rna_type, data):
    """
    This can ignore any conflict where the sequence probably comes from a
    mitochondria but it matches a bacterial rRNA. In that case we do not
    warn since this is expected from evolution.
    """
    return data.is_mitochondrial() and \
        'rRNA' == rna_type and \
        data.rfam_hits[0].model in {
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

    hit = data.rfam_hits[0]
    if not hit.model_domain:
        return False

    if not data.domains():
        return False

    return hit.model_domain not in data.domains() and \
        not is_ignorable_mito_conflict(rna_type, data)


def incomplete_sequence(_, data):
    """
    Detect if the given sequence is incomplete. This will work by checking if
    the Rfam hits cover enough of the sequence and the hit is complete enough.
    """

    hits = data.rfam_hits
    if len(hits) != 1:
        return False

    hit = hits[0]
    if hit.model not in INCOMPLETE_FAMILIES:
        return False

    return hit.model_info.completeness <= 0.5 and \
        hit.sequence_info.completeness >= 0.9


def missing_rfam_match(rna_type, data):
    """
    Detect if the given sequence, with the given RNA type is missing an Rfam
    match.
    """

    if rna_type not in EXPECTED_MATCHES:
        return False

    if rna_type == 'tRNA' and \
            data.is_mitochondrial() and \
            not data.rfam_hits:
        return False

    required = EXPECTED_MATCHES[rna_type]
    hits = {h.model for h in data.rfam_hits}
    return not hits.intersection(required)


def mismatching_rna_type(rna_type, data):
    """
    This will raise a warning if the hits have different RNA types than the
    computed type. This indicates some annotation issue. Currently this only
    works for miRNA's.
    """

    if rna_type != 'miRNA':
        return None
    rna_types = {hit.model_rna_type for hit in data.rfam_hits}
    return rna_types != {rna_type}
