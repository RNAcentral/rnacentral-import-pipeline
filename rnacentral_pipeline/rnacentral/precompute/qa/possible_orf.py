# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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


NAME = "possible_orf"


def validate(sequence: Sequence) -> QaResult:
    if not sequence.orf_info:
        return QaResult.ok(NAME)

    sources = sequence.orf_info.all_sources()
    count = len(sources)
    sources = ", ".join(sources)
    message = f"This sequence contains a possible orf, as annotated by {sources}"
    if count > 1:
        message = f"This sequence contains {count} possible orfs, as annotated by {sources}"
    return QaResult.not_ok(NAME, message)
