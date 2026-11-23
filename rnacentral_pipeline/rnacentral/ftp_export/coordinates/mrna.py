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

import csv
import itertools as it
import json
import sys
import typing as ty
from collections import OrderedDict

from gffutils import Feature

from .gff3 import write_gff_text


def within(exons, start: int, stop: int) -> ty.List[ty.Tuple[int, int]]:
    return [(max(s, start), min(e, stop)) for s, e in exons if s <= stop and e >= start]


def features(record: dict) -> ty.Iterable[Feature]:
    """
    Turn one mRNA from mrna.sql into a GFF3 gene model: the mRNA, its exons,
    and the UTR and CDS parts of those exons, split at the stored CDS span.
    """
    exons = [(e["start"], e["stop"]) for e in record["exons"]]
    note = json.loads(record["note"])
    cds_start, cds_stop = note["cds"]
    strand = "+" if record["strand"] > 0 else "-"
    transcript = record["transcript"]

    def feature(kind, start, stop, attributes, frame="."):
        return Feature(
            seqid=record["chromosome"],
            source="RNAcentral",
            featuretype=kind,
            start=start,
            end=stop,
            strand=strand,
            frame=frame,
            attributes=attributes,
        )

    attributes = OrderedDict(
        [
            ("ID", [transcript]),
            ("Name", [record["rna_id"]]),
            ("transcript_id", [transcript]),
            ("gene_id", [record["gene"]]),
        ]
    )
    if record["gene_name"]:
        attributes["gene_name"] = [record["gene_name"]]
    if note.get("tags"):
        attributes["tag"] = note["tags"]
    yield feature("mRNA", exons[0][0], exons[-1][1], attributes)

    parent = OrderedDict([("Parent", [transcript])])
    for start, stop in exons:
        yield feature("exon", start, stop, parent)

    upstream = within(exons, exons[0][0], cds_start - 1)
    downstream = within(exons, cds_stop + 1, exons[-1][1])
    five, three = (upstream, downstream) if strand == "+" else (downstream, upstream)

    def utr(key):
        attributes = OrderedDict(parent)
        if record[key]:
            attributes["Name"] = [record[key]]
        return attributes

    for start, stop in five:
        yield feature("five_prime_UTR", start, stop, utr("five_prime_utr"))

    # A CDS with an incomplete start begins cds_phase bases into a codon
    first_phase = note.get("cds_phase", 0)
    cds = within(exons, cds_start, cds_stop)
    done = 0
    for start, stop in cds if strand == "+" else reversed(cds):
        yield feature("CDS", start, stop, parent, frame=str((first_phase - done) % 3))
        done += stop - start + 1

    for start, stop in three:
        yield feature("three_prime_UTR", start, stop, utr("three_prime_utr"))


def from_file(handle, output, allow_none=False):
    csv.field_size_limit(sys.maxsize)
    records = (json.loads(row[0]) for row in csv.reader(handle))
    found = it.chain.from_iterable(features(r) for r in records)
    write_gff_text(found, output, allow_no_features=allow_none)
