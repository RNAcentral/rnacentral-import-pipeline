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
import typing as ty

from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.silva import helpers


def parse(handle, taxonomy_path) -> ty.Iterable[Entry]:
    taxonomy = SqliteDict(filename=taxonomy_path)
    reader = csv.DictReader(handle, delimiter="\t")
    entries = map(lambda r: helpers.as_entry(taxonomy, r), reader)
    return filter(None, entries)
