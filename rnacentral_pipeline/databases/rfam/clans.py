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
from pathlib import Path

from rnacentral_pipeline import schemas
from rnacentral_pipeline.parquet_writers import row_writer

# Rfam's TSV column names, in the order load_rfam_clans expects them.
FIELDS = [
    "id",
    "name",
    "description",
    "family_count",
]


def from_file(filename, output):
    reader = csv.DictReader(filename, delimiter="\t")
    rows = ([row[f] for f in FIELDS] for row in reader)
    with row_writer(Path(output), schemas.RFAM_CLANS) as writer:
        writer.writerows(rows)
