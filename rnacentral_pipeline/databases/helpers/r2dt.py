# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

import re
import typing as ty
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Both CRW and RiboVision ship bpseq2fasta output, whose header names the bpseq
# it came from. The paths differ between the two, so match only the basename.
MODEL_ID = re.compile(r"([^/\s]+)\.bpseq")


def fasta_entries(directories: ty.Iterable[Path]) -> ty.Iterable[SeqRecord]:
    """
    Read the model sequences R2DT ships as bpseq2fasta output: a header naming
    the source bpseq, the sequence, then the dot-bracket structure. The
    structure is dropped -- only the sequence is imported.
    """
    for directory in directories:
        for path in sorted(directory.glob("*.fasta")):
            with path.open("r") as raw:
                lines = [line.strip() for line in raw if line.strip()]

            if len(lines) < 2:
                raise ValueError(f"{path} is not bpseq2fasta output")

            match = MODEL_ID.search(lines[0])
            if match is None:
                raise ValueError(f"Could not get model id from {lines[0]}")

            yield SeqRecord(Seq(lines[1]), id=match.group(1))
