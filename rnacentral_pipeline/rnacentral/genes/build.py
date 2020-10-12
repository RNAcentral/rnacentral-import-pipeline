# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
import operator as op

from rnacentral_pipeline import psql
from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import \
    write_bed_text

from . import data, rrna


def load(handle):
    for entry in psql.json_handler(handle):
        yield data.UnboundLocation.build(entry)


def split_key(value: data.UnboundLocation):
    extent = value.extent
    return (extent.chromosome, extent.strand, extent.insdc_rna_type)


def split(handle):
    entries = load(data)
    for (key, locations) in it.groupby(entries, key):
        if key[-1] != "rRNA":
            continue

        yield list(locations)


def build(data):
    for group in data:
        genes = rrna.build(group)
        for gene in genes:
            for writeable in gene.writeable():
                yield writeable


def write_genes(raw, output):
    data = split(raw)
    writer = csv.writer(output)
    writer.writerows(build(data))


def write_bed(raw, output):
    data = split(raw)
    bed = map(d.as_bed() for d in data)
    write_bed_text(bed, output)


def write(raw, output, as_bed):
    if as_bed:
        write_bed(raw, output)
    else:
        write_genes(raw, output)
