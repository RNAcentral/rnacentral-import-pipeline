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

ALLOWED_MITO_FAMILIES = {
    'RF00177',  # Bacterial small subunit ribosomal RNA
    'RF02541',  # Bacterial large subunit ribosomal RNA
    'RF01959',  # Archaeal small subunit ribosomal RNA
    'RF02540',  # Archaeal large subunit ribosomal RNA
}


def is_ignorable_mito_conflict(rna_type, data):
    """
    This can ignore any conflict where the sequence probably comes from a
    mitochondria but it matches a bacterial rRNA. In that case we do not
    warn since this is expected from evolution.
    """
    return data.is_mitochondrial() and \
        'rRNA' == rna_type and \
        data.rfam_hits[0].model in ALLOWED_MITO_FAMILIES


class Validator(object):
    """
    Validator for checking for possible sequence contamination.
    """
    name = 'possible_contamination'

    def status(self, rna_type, data):
        """
        Check if the given sequence with the given RNA type is likely
        contamination. This checks if the domain of the sequence and the domain
        of the hits disagree.
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

    def message(self, _, data):
        """
        Produce a warning message about the issue detected. This assumes that
        there was a warning.
        """

        return (
            'This {sequence_name} sequence matches a {match_domain} '
            'Rfam model: <a href="{url}">{model_name}</a>. '
        ).format(
            sequence_name=data,
            match_domain=None,
            model_name=None,
        )
