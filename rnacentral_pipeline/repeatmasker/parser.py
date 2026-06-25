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

import typing as ty
from collections import defaultdict
from pathlib import Path

from rnacentral_pipeline.repeatmasker.data import (
    RepeatmaskerResult,
    RepeatRegion,
)
from rnacentral_pipeline.stopfree.scan import read_fasta_records


def parse_out(handle: ty.TextIO) -> ty.Iterator[RepeatRegion]:
    """
    Parse the hit lines of a RepeatMasker .out file. The file has a 2-line
    header plus a blank line; data lines are whitespace separated and start with
    the integer Smith-Waterman score. Columns: score, %div, %del, %ins, query,
    begin, end, (left), strand, repeat, class/family, ...
    """
    for raw in handle:
        fields = raw.split()
        if len(fields) < 11 or not fields[0].isdigit():
            continue

        strand = "-" if fields[8] == "C" else "+"
        yield RepeatRegion(
            urs_taxid=fields[4],
            start=int(fields[5]) - 1,
            stop=int(fields[6]),
            repeat_name=fields[9],
            repeat_class=fields[10],
            strand=strand,
            divergence=float(fields[1]),
            score=int(fields[0]),
        )


def parse(
    sequences: Path, repeatmasker_out: Path
) -> ty.Iterator[ty.Tuple[RepeatmaskerResult, ty.List[RepeatRegion]]]:
    """
    Build per-URS_taxid results from the scanned FASTA and the RepeatMasker .out
    file. A result row is emitted for every input sequence (including those with
    no repeats) so the tracking table records absence and the workflow stays
    idempotent.
    """
    regions: ty.Dict[str, ty.List[RepeatRegion]] = defaultdict(list)
    with repeatmasker_out.open("r") as handle:
        for region in parse_out(handle):
            regions[region.urs_taxid].append(region)

    for urs_taxid, sequence in read_fasta_records(sequences):
        hits = regions.get(urs_taxid, [])
        result = RepeatmaskerResult.build(urs_taxid, len(sequence), hits)
        yield result, hits
