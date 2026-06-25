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

from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult

NAME = "from_repetitive_region"


def validate(sequence: Sequence) -> QaResult:
    if sequence.from_repeat is None:
        return QaResult.null(NAME)
    if not sequence.from_repeat:
        return QaResult.ok(NAME)
    message = "This sequence overlaps a repetitive region, as annotated by RepeatMasker"
    return QaResult.not_ok(NAME, message)
