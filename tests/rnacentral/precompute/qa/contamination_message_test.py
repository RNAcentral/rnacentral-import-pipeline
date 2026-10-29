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

from rnacentral_pipeline.rnacentral.precompute.qa import contamination

from .. import builders as b

RRNA = "SO:0000252"


def test_message_does_not_crash_with_no_common_name_or_species():
    """
    contamination.message() only assigns its local `sequence_name` inside
    the `len(common_name) == 1` branch or the nested `if species:` branch.
    With accessions that give neither a single common_name nor any species,
    neither branch runs and `sequence_name` is referenced unbound. This
    covers exactly that: no accessions at all, so both `common_name` and
    `species` come out empty.
    """
    sequence = b.sequence(
        accessions=[],
        rfam_hits=[b.rfam_hit("RF00001", RRNA, model_domain="Bacteria")],
    )
    message = contamination.message(sequence)
    assert isinstance(message, str)
