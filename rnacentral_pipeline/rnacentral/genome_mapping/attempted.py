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

from rnacentral_pipeline import psql


def parse(handle, taxid, assembly_id):
    for entry in psql.json_handler(handle):
        urs_taxid = f"{entry['id']}_{taxid}"
        yield [urs_taxid, assembly_id]


def write(handle, taxid, assembly_id, output):
    writer = csv.writer(output)
    writer.writerows(parse(handle, assembly_id))
