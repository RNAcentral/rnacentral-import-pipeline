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

import json


def split(handle):
    output = Path(out)
    key = op.itemgetter("rna_type", "chromosome", "strand")
    entries = sorted(psql.json_handle(handle), key=key)
    for (key, locations) in it.groupby(entries, key):
        if key[0] != "rRNA":
            continue

        yield list(locations)


def build(data):
    for group in data:
        for entry in group:
            yield [
                entry["taxid"],
                entry["assembly_id"],
                entry["chromosome"],
                entry["strand"],
                entry["exons"][0]["exon_start"],
                entry["exons"][-1]["exon_stop"],
                entry["urs_taxid"],
                entry["region_id"],
                entry["qa"]["has_issue"],
            ]


def write_genes(raw, output):
    data = split(raw)
    writer = csv.writer(output)
    writer.writerows(build(data))
