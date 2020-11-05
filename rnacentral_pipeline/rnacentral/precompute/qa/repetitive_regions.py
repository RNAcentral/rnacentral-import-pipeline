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
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult


def validate(context: Context, rna_type: str, sequence: Sequence) -> QaResult:
    repeats = context.repeats
    for coordinate in sequence.coordinates:
        if not repeats.has_assembly(coordinate.assembly_id):
            continue
        enclosed = repeats.envelops(
            coordinate.assembly_id,
            coordinate.chromosome,
            coordinate.start,
            coordinate.stop,
        )
        if enclosed:
            return QaResult.not_ok(
                "from_repetitive_region", "This sequence overlaps a repetitive region"
            )
    return QaResult.ok("from_repetitive_region")
