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

from rnacentral_pipeline.tcode.data import TcodeResult


def parse(tcode_output: Path) -> ty.Iterable[TcodeResult]:
    sequence = None
    size: ty.Optional[int] = None
    scores: ty.List[float] = []

    def flush() -> ty.Optional[TcodeResult]:
        if not sequence:
            return None
        if scores:
            mean_score = statistics.mean(scores)
            std_score = statistics.stdev(scores) if len(scores) > 1 else 0.0
        else:
            mean_score = None
            std_score = None
        return TcodeResult.build(sequence, size, mean_score, std_score)

    with tcode_output.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("# Sequence:"):
                result = flush()
                if result:
                    yield result
                scores = []
                size = None
                sequence = line.split(":", 1)[1].strip().split()[0]
                continue

            if line.startswith("# Total_length:"):
                try:
                    size = int(line.split(":", 1)[1].strip())
                except ValueError:
                    size = None
                continue

            if line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 4:
                continue

            try:
                score = float(fields[3])
            except ValueError:
                continue

            scores.append(score)

    result = flush()
    if result:
        yield result
