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


def _split_taxid(seq_id: str) -> ty.Tuple[str, bool]:
    if "_" not in seq_id:
        return seq_id, False
    base, suffix = seq_id.rsplit("_", 1)
    return (base, True) if suffix.isdigit() else (seq_id, False)


def _build_result(
    sequence: str,
    size: ty.Optional[int],
    header_len: ty.Optional[int],
    scores: ty.List[float],
) -> TcodeResult:
    if scores:
        mean_score = statistics.mean(scores)
        std_score = statistics.stdev(scores) if len(scores) > 1 else 0.0
    else:
        mean_score = None
        std_score = None
    final_size = size if size is not None else header_len
    return TcodeResult.build(sequence, final_size, mean_score, std_score)


def parse(tcode_output: Path) -> ty.Iterable[TcodeResult]:
    sequence = None
    size: ty.Optional[int] = None
    header_len: ty.Optional[int] = None
    scores: ty.List[float] = []
    results_by_base: ty.Dict[str, TcodeResult] = {}
    seq_pattern = re.compile(r"^#?\s*Sequence:\s+(\S+)\s+from:\s+(\d+)\s+to:\s+(\d+)")

    with tcode_output.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            match = seq_pattern.match(line)
            if match:
                if sequence:
                    result = _build_result(sequence, size, header_len, scores)
                    base, has_taxid = _split_taxid(result.urs)
                    if has_taxid:
                        results_by_base[base] = result
                scores = []
                size = None
                header_len = None
                sequence = match.group(1)
                try:
                    header_len = int(match.group(3)) - int(match.group(2)) + 1
                except ValueError:
                    header_len = None
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
                scores.append(float(fields[3]))
            except ValueError:
                continue

    if sequence:
        result = _build_result(sequence, size, header_len, scores)
        base, has_taxid = _split_taxid(result.urs)
        if has_taxid:
            results_by_base[base] = result

    for result in results_by_base.values():
        yield result
