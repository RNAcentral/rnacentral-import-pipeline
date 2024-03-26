# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This contains the logic for parsing MGI data files and producing Entry objects
for export to usable flat files.
"""

import csv
import itertools as it
import operator as op

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.writers import build_entry_writer

from . import helpers


def lines(raw):
    """
    Produces an iterable of all ines in the file. This will correct the issues
    with header being over 2 lines so a normal CSV parser can parse the file.
    """

    header = "\t".join([next(raw).strip(), next(raw).strip()])
    header = header.lower()
    yield header.replace(" ", "_")
    for line in raw:
        yield line


def as_entry(data):
    yield Entry(
        primary_id=helpers.primary_id(data),
        accession=helpers.accession(data),
        ncbi_tax_id=helpers.taxon_id(data),
        database="MGI",
        sequence="",
        exons=[],
        rna_type=helpers.infer_rna_type(data) or "",
        url="",
        xref_data=helpers.xref_data(data),
        chromosome=helpers.chromosome(data),
        species=helpers.species(data),
        common_name=helpers.common_name(data),
        lineage=helpers.lineage(data),
        gene=helpers.gene(data),
        optional_id=helpers.symbol(data),
        description=helpers.name(data),
        seq_version="1",
        location_start=helpers.start(data),
        location_end=helpers.stop(data),
        references=helpers.references(data),
    )


def parse(raw):
    """
    Parses the file and produces an iterable of all entries as MGI objects.
    """

    data = csv.DictReader(lines(raw), delimiter="\t")
    data = map(as_entry, data)
    return filter(op.attrgetter("rna_type"), data)


def from_file(raw, output):
    writer = build_entry_writer(parse)
    writer(output, raw)
