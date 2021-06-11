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
from rnacentral_pipeline.rnacentral.precompute.qa.data import QaResult

EXPECTED_MATCHES = {
    "rRNA": {
        "RF00001",  # 5S ribosomal RNA
        "RF00002",  # 5.8S ribosomal RNA
        "RF00177",  # Bacterial small subunit ribosomal RNA
        "RF01959",  # Archaeal small subunit ribosomal RNA
        "RF01960",  # Eukaryotic small subunit ribosomal RNA
        "RF02540",  # Archaeal large subunit ribosomal RNA
        "RF02541",  # Bacterial large subunit ribosomal RNA
        "RF02542",  # Microsporidia small subunit ribosomal RNA
        "RF02543",  # Eukaryotic large subunit ribosomal RNA
        "RF02547",  # mito 5S RNA
    },
    "tRNA": {
        "RF00005",
        "RF01852",
    },  # tRNA  # Selenocysteine tRNA
}


def href(model_id):
    return f'<a href="http://rfam.org/family/{model_id}">{model_id}</a>'


def validate(context: Context, rna_type: str, sequence: Sequence) -> QaResult:
    if rna_type not in EXPECTED_MATCHES:
        return QaResult.ok("missing_rfam_match")

    if rna_type == "tRNA" and sequence.is_mitochondrial() and not sequence.rfam_hits:
        return QaResult.ok("missing_rfam_match")

    required = EXPECTED_MATCHES[rna_type]
    hits = {h.model for h in sequence.rfam_hits}
    if not hits.intersection(required):
        possible = sorted(EXPECTED_MATCHES[rna_type])

        article = "the"
        if len(possible) > 1:
            article = "a"

        models = [href(p) for p in sorted(possible)]
        expected = ", ".join(models)
        message = f"No match to {article} {rna_type} Rfam model ({expected})"
        return QaResult.not_ok("missing_rfam_match", message)
    return QaResult.ok("missing_rfam_match")
