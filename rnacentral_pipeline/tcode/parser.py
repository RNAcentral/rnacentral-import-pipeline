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

from pathlib import Path
import statistics
import typing as ty
import re

from rnacentral_pipeline.tcode.data import TcodeResult


def _build_result(
    sequence: str,
    length: ty.Optional[int],
    scores: ty.List[float],
) -> TcodeResult:
    if scores:
        mean_score = statistics.mean(scores)
        std_score = statistics.stdev(scores) if len(scores) > 1 else 0.0
    else:
        mean_score = None
        std_score = None
    return TcodeResult.build(sequence, length, mean_score, std_score)


def parse(tcode_output: Path) -> ty.Iterable[TcodeResult]:
    """
    Assumes input sequences are pre-filtered ( >200 bp to avoid tcode crashing or
    unreliable results). 
    For each sequence block, parse scores and emit mean/std. If there is one score,
    mean is that score and std is 0.0. If no scores are reported, mean/std are NaN.
    """
    sequence = None
    length: ty.Optional[int] = None
    scores: ty.List[float] = []
    seq_pattern = re.compile(r"^#?\s*Sequence:\s+(\S+)\s+from:\s+(\d+)\s+to:\s+(\d+)")

    with tcode_output.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            match = seq_pattern.match(line)
            if match:
                if sequence:
                    yield _build_result(sequence, length, scores)
                scores = []
                try:
                    length = int(match.group(3)) - int(match.group(2)) + 1
                except ValueError:
                    length = None
                sequence = match.group(1)
                continue

            if line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 4:
                continue

            try:
                scores.append(float(fields[3]))
            except ValueError:
                continue

    if sequence:
        yield _build_result(sequence, length, scores)
