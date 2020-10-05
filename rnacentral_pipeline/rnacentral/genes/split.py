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

import itertools as it
import operator as op

from rnacentral_pipeline import psql


def write_split(handle, out):
    output = Path(out)
    key = op.itemgetter("rna_type", "chromosome")
    entries = sorted(psql.json_handle(handle), key=key)
    for ((rna_type, chromosome), locations) in it.groupby(entries, key):
        if rna_type != "rRNA":
            continue

        path = output / f"{rna_type}-{chromosome}.json"
        with path.open("w") as out:
            json.dump(list(locations), out)
