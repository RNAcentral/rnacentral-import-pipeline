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
import logging
import operator as op

from Bio import SeqIO

from rnacentral_pipeline import psql

LOGGER = logging.getLogger(__name__)


id_attribute = op.attrgetter("id")
id_key = op.itemgetter("id")


def append_taxid(taxid, getter=id_key):
    def fn(entry):
        urs = getter(entry)
        return f"{urs}_{taxid}"

    return fn


def json_parser(handle, id_generator=id_key, extra_fields=[]):
    for entry in psql.json_handler(handle):
        yield [id_generator(entry)] + extra_fields


def fasta_parser(handle, id_generator=id_attribute, extra_fields=[]):
    for record in SeqIO.parse(handle, "fasta"):
        seq_id = id_generator(record)
        if not seq_id:
            raise ValueError("Failed to find seq id for %s" % record)
        yield [seq_id] + extra_fields


def parse_rfam_version(handle):
    for line in handle:
        if line.startswith("Release"):
            parts = line.split(" ")
            return parts[1].strip()
    raise ValueError("Could not find version in file")


def write(data, output, require_attempt=True):
    writer = csv.writer(output)
    seen = False
    for row in data:
        writer.writerows(row)
        seen = True
    if not seen:
        LOGGER.error("Nothing was attempted")
        if require_attempt:
            raise ValueError("Found nothing was attempted")


def genome_mapping(handle, assembly_id, output):
    data = fasta_parser(handle, extra_fields=[assembly_id])
    write(data, output)


def qa(handle, name, version_file, output):
    if name == "rfam":
        version = parse_rfam_version(version_file)
    else:
        raise ValueError(f"Unknown QA type: {name}")
    data = fasta_parser(handle, extra_fields=[name, version])
    write(data, output)


def r2dt(handle, output):
    data = fasta_parser(handle)
    write(data, output)
