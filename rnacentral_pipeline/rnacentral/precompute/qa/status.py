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


from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.qa import contamination
from rnacentral_pipeline.rnacentral.precompute.qa import (
    incomplete_sequence as incomplete,
)
from rnacentral_pipeline.rnacentral.precompute.qa import missing_rfam_match as missing
from rnacentral_pipeline.rnacentral.precompute.qa import (
    repetitive_regions as repetitive,
)
from rnacentral_pipeline.rnacentral.precompute.qa import possible_orf
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaStatus


def status(context: Context, sequence: Sequence, rna_type: str) -> QaStatus:
    """
    Generate the QaStatus for a given Sequence object
    """

    return QaStatus(
        incomplete_sequence=incomplete.validate(sequence),
        possible_contamination=contamination.validate(rna_type, sequence),
        missing_rfam_match=missing.validate(rna_type, sequence),
        from_repetitive_region=repetitive.validate(context, sequence),
        possible_orf=possible_orf.validate(sequence),
    )
