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

import typing as ty

from rnacentral_pipeline.rnacentral.precompute.data.orf import OrfInfo
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult

# These must stay exactly these strings - they are QaResult.name and match
# the qa.ctl load columns / QaStatus fields.
NAME_CPAT = "possible_orf"
NAME_STOPFREE = "possible_orf_stopfree"
NAME_TCODE = "possible_orf_tcode"


def _annotated_by(sources: ty.List[str]) -> str:
    count = len(sources)
    joined = ", ".join(sources)
    if count > 1:
        return f"This sequence contains {count} possible ORFs, as annotated by {joined}"
    return f"This sequence contains a possible ORF, as annotated by {joined}"


def _cpat_message(orf_info: ty.Optional[OrfInfo]) -> str:
    # ponytail: CPAT can report actual ORF coordinates, unlike stopfree/tcode
    # which are bare booleans - this is where a coordinate string would get
    # appended once OrfInfo actually carries one. Nothing populates that today.
    if orf_info is None:
        return _annotated_by(["CPAT"])
    return _annotated_by(orf_info.all_sources())


_STOPFREE_MESSAGE = _annotated_by(["stopfree"])
_TCODE_MESSAGE = _annotated_by(["tcode"])


def _result(name: str, flag: ty.Optional[bool], message: str) -> QaResult:
    if flag is None:
        return QaResult.null(name)
    if not flag:
        return QaResult.ok(name)
    return QaResult.not_ok(name, message)


def validate(sequence: Sequence) -> ty.Tuple[QaResult, QaResult, QaResult]:
    """
    Check all three coding-potential signals for a sequence. Returns
    (cpat, stopfree, tcode), matching the three fixed QaStatus fields /
    qa.ctl columns they get assigned to.
    """

    cpat = _result(NAME_CPAT, sequence.possible_orf, _cpat_message(sequence.orf_info))
    stopfree = _result(NAME_STOPFREE, sequence.possible_orf_stopfree, _STOPFREE_MESSAGE)
    tcode = _result(NAME_TCODE, sequence.possible_orf_tcode, _TCODE_MESSAGE)
    return (cpat, stopfree, tcode)
