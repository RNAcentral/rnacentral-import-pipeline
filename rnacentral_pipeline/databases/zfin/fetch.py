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

import gzip
import json
import urllib.request as request
from contextlib import closing


def fetch(url):
    with closing(request.urlopen(url)) as compressed:
        with gzip.GzipFile(None, 'rb', 9, compressed) as raw:
            data = json.load(raw)

    # Fix weird PMID formatting
    pubs = []
    for publication in data['metaData']['publications']:
        pubs.append(publication.replace(': ', ':'))
    data['metaData']['publications'] = pubs

    # Strip out entries without SO terms:
    valid = []
    for entry in data['data']:
        if 'soTermId' not in entry:
            continue
        valid.append(entry)
    data['data'] = valid
    return data
