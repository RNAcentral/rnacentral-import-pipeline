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

from rnacentral_pipeline.databases.data import CODING_DATABASES
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.qa import coding_potential, contamination
from rnacentral_pipeline.rnacentral.precompute.qa import (
    incomplete_sequence as incomplete,
)
from rnacentral_pipeline.rnacentral.precompute.qa import missing_rfam_match as missing
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult, QaStatus


def status(context: Context, sequence: Sequence, rna_type: str) -> QaStatus:
    """
    Generate the QaStatus for a given Sequence object
    """

    # These checks judge ncRNAs; on mRNAs they only raise meaningless flags, so a
    # sequence known only as a coding database's mRNA is not checked.
    accessions = sequence.accessions
    if accessions and all(
        a.database in CODING_DATABASES and a.rna_type.so_term.so_id == "SO:0000234"
        for a in accessions
    ):
        # The QA table needs the first four checks set; only the ORF columns take nulls.
        ok, null = QaResult.ok, QaResult.null
        return QaStatus(
            incomplete_sequence=ok("incomplete_sequence"),
            possible_contamination=ok("possible_contamination"),
            missing_rfam_match=ok("missing_rfam_match"),
            from_repetitive_region=ok("from_repetitive_region"),
            possible_orf=null("possible_orf"),
            possible_orf_stopfree=null("possible_orf_stopfree"),
            possible_orf_tcode=null("possible_orf_tcode"),
        )

    orf, stopfree, tcode = coding_potential.validate(sequence)
    return QaStatus(
        incomplete_sequence=incomplete.validate(sequence),
        possible_contamination=contamination.validate(rna_type, sequence),
        missing_rfam_match=missing.validate(rna_type, sequence),
        # ponytail: repeat check removed (it never had data, see
        # docs/repetitive-region-qa-disabled.md); the column stays until RepeatMasker
        # QA replaces it.
        from_repetitive_region=QaResult.ok("from_repetitive_region"),
        possible_orf=orf,
        possible_orf_stopfree=stopfree,
        possible_orf_tcode=tcode,
    )
