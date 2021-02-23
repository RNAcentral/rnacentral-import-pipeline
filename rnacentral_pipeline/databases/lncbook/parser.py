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
from io import StringIO

from rnacentral_pipeline.databases.generic import parser as generic


def parse(handle):
    raw = json.load(handle)
    data = []
    for ncrna in raw['data']:
        regions = ncrna['genomeLocactions']
        regions = filter(lambda r: r['assembly'] == 'GRCh38', regions)
        regions = list(regions)
        if not regions:
            continue
        ncrna['genomeLocactions'] = regions
        data.append(ncrna)
    if not data:
        raise ValueError("All ncRNA are not from GRCh38, failing")
    raw = {
        "metaData": raw["metaData"],
        "data": data,
    }
    processed = StringIO()
    json.dump(raw, processed)
    processed.seek(0)

    for entry in  generic.parse(processed):
        updates = []
        for region in entry.regions:
            if region.chromosome == 'M':
                updates.append(attr.evolve(region, chromosome='MT'))
            else:
                updates.append(region)
        entry = attr.evolve(entry, regions=updates)
        yield entry
