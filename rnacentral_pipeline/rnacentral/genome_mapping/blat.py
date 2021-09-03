from __future__ import annotations

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

import csv
import operator as op
import itertools as it
import logging
import typing as ty

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline import utils
from rnacentral_pipeline.databases.data.regions import Exon
from rnacentral_pipeline.databases.data.regions import SequenceRegion
from rnacentral_pipeline.databases.data.regions import CoordinateSystem

LOGGER = logging.getLogger(__name__)

FIELDS = [
    "matches",  # Number of bases that match that aren't repeats
    "misMatches",  # Number of bases that don't match
    "repMatches",  # Number of bases that match but are part of repeats
    "nCount",  # Number of "N" bases
    "qNumInsert",  # Number of inserts in query
    "qBaseInsert",  # Number of bases inserted in query
    "tNumInsert",  # Number of inserts in target
    "tBaseInsert",  # Number of bases inserted in target
    "strand",  # "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
    "qName",  # Query sequence name
    "qSize",  # Query sequence size.
    "qStart",  # Alignment start position in query
    "qEnd",  # Alignment end position in query
    "tName",  # Target sequence name
    "tSize",  # Target sequence size
    "tStart",  # Alignment start position in target
    "tEnd",  # Alignment end position in target
    "blockCount",  # Number of blocks in the alignment (a block contains no gaps)
    "blockSizes",  # Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
    "qStarts",  # Comma-separated list of starting positions of each block in query
    "tStarts",  # Comma-separated list of starting positions of each block in target
]


@attr.s(frozen=True)
class BlatHit(object):
    upi = attr.ib(validator=is_a(str), converter=str)
    sequence_length = attr.ib(validator=is_a(int))
    matches = attr.ib(validator=is_a(int))
    target_insertions = attr.ib(validator=is_a(int))
    region = attr.ib(validator=is_a(SequenceRegion))

    @classmethod
    def build(cls, assembly_id: str, raw: ty.Dict[str, ty.Any]) -> BlatHit:
        parts = zip(raw["tStarts"], raw["blockSizes"])
        exons = [Exon(s, s + l) for (s, l) in parts]

        return cls(
            upi=raw["qName"],
            sequence_length=raw["qSize"],
            matches=raw["matches"],
            target_insertions=raw["tBaseInsert"],
            region=SequenceRegion(
                assembly_id=assembly_id,
                chromosome=raw["tName"],
                strand=raw["strand"],
                exons=exons,
                coordinate_system=CoordinateSystem.zero_based(),
            ),
        )

    @property
    def name(self) -> str:
        return self.region.name(upi=self.upi)

    @property
    def match_fraction(self) -> float:
        return float(self.matches) / float(self.sequence_length)

    def writeable(self):
        return self.region.writeable(self.upi, is_upi=True)


def select_possible(hit: BlatHit) -> bool:
    if hit.matches < 100 and hit.target_insertions > 25:
        return False
    if hit.matches == hit.sequence_length:
        return True
    if (
        hit.sequence_length > 15
        and hit.match_fraction > 0.95
        and hit.match_fraction < 1
    ):
        return True
    return False


def select_best(hits: ty.Iterable[BlatHit]) -> ty.List[BlatHit]:
    hits = list(hits)
    best = max(hits, key=op.attrgetter("match_fraction"))
    return [h for h in hits if h.match_fraction >= best.match_fraction]


def parse_psl(assembly_id: str, handle: ty.IO) -> ty.Iterable[BlatHit]:
    to_split = ["blockSizes", "qStarts", "tStarts"]
    for row in csv.reader(handle, delimiter="\t"):
        result: ty.Dict[str, ty.Any] = dict(zip(FIELDS, row))
        for key in to_split:
            result[key] = [int(v) for v in result[key].split(",") if v]
        lens = {len(result[k]) for k in to_split}
        assert len(lens) == 1

        for key, value in result.items():
            if key not in to_split and "Name" not in key and key != "strand":
                result[key] = int(value)

        yield BlatHit.build(assembly_id, result)


def select_hits(hits: ty.Iterable[BlatHit], sort=False) -> ty.Iterable[BlatHit]:
    key = op.attrgetter("upi")
    if sort:
        hits = sorted(hits, key=key)

    for upi, subhits in it.groupby(hits, key=key):
        selected = list(filter(select_possible, subhits))

        if not selected:
            LOGGER.warn("No possible matches for %s", upi)
            continue

        best = select_best(selected)
        if not best:
            raise ValueError("Failed to select a best hit for %s" % upi)

        for hit in best:
            yield hit


def write_importable(handle, output):
    hits = utils.unpickle_stream(handle)
    writeable = map(op.methodcaller("writeable"), hits)
    writeable = it.chain.from_iterable(writeable)
    csv.writer(output).writerows(writeable)


def as_pickle(assembly_id, hits, output):
    parsed = parse_psl(assembly_id, hits)
    utils.pickle_stream(parsed, output)


def select_pickle(handle, output, sort=False):
    hits = utils.unpickle_stream(handle)
    selected = select_hits(hits, sort=sort)
    utils.pickle_stream(selected, output)
