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

from pathlib import Path
import typing as ty

from Bio import SeqIO

from rnacentral_pipeline import psql

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.crw import helpers


def parse(metadata_handle: ty.IO, sequences: Path) -> ty.Iterable[Entry]:
    indexed = SeqIO.index(sequences, "fasta")
    for entry in psql.json_handler(metadata_handle):
        data = helpers.as_entry(entry, indexed)
        if data:
            yield data
