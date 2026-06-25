from __future__ import annotations

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

import json
import typing as ty

import attr
from attr.validators import instance_of as is_a


@attr.s()
class RepeatRegion:
    """A single repeat hit from a RepeatMasker .out file."""

    urs_taxid = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int))  # 0-based, inclusive
    stop = attr.ib(validator=is_a(int))  # inclusive
    repeat_name = attr.ib(validator=is_a(str))
    repeat_class = attr.ib(validator=is_a(str))
    strand = attr.ib(validator=is_a(str))  # '+' or '-'
    divergence = attr.ib(validator=is_a(float))
    score = attr.ib(validator=is_a(int))

    @property
    def length(self) -> int:
        return self.stop - self.start

    def writeable(self) -> ty.List[str]:
        urs, taxid = self.urs_taxid.split("_", 1)
        metadata = {
            "repeat": self.repeat_name,
            "class": self.repeat_class,
            "strand": self.strand,
            "divergence": self.divergence,
            "score": self.score,
        }
        return [urs, taxid, str(self.start), str(self.stop), json.dumps(metadata)]


@attr.s()
class RepeatmaskerResult:
    """Per-URS_taxid presence summary, mirroring tcode/stopfree/cpat results."""

    urs_taxid = attr.ib(validator=is_a(str))
    has_repeats = attr.ib(validator=is_a(bool))
    repeat_coverage = attr.ib(validator=is_a(float))
    repeat_count = attr.ib(validator=is_a(int))

    @classmethod
    def build(
        cls, urs_taxid: str, length: int, regions: ty.List[RepeatRegion]
    ) -> RepeatmaskerResult:
        # Naive covered-base count; repeat regions from RepeatMasker do not
        # overlap within a query, so summing hit lengths is safe enough here.
        # ponytail: naive sum, switch to interval merge if overlaps ever appear.
        covered = sum(r.length for r in regions)
        coverage = covered / length if length else 0.0
        return cls(
            urs_taxid=urs_taxid,
            has_repeats=bool(regions),
            repeat_coverage=round(coverage, 4),
            repeat_count=len(regions),
        )

    def writeable(self) -> ty.List[str]:
        return [
            self.urs_taxid,
            str(self.has_repeats),
            str(self.repeat_coverage),
            str(self.repeat_count),
        ]


@attr.s()
class RepeatmaskerWriter:
    results = attr.ib()
    features = attr.ib()

    def write(
        self,
        data: ty.Iterable[ty.Tuple[RepeatmaskerResult, ty.List[RepeatRegion]]],
    ):
        for result, regions in data:
            self.results.writerow(result.writeable())
            for region in regions:
                self.features.writerow(region.writeable())
