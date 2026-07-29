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
from time import sleep
import urllib.error as error
import urllib.request as request
from contextlib import closing

TRANSIENT_HTTP_STATUS_CODES = {429, 500, 502, 503, 504}
MAX_RETRIES = 5


def fetch(url):
    data = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with closing(request.urlopen(url, timeout=30)) as raw:
                data = json.load(raw)
            break
        except error.HTTPError as err:
            if attempt == MAX_RETRIES or err.code not in TRANSIENT_HTTP_STATUS_CODES:
                raise
        except error.URLError:
            if attempt == MAX_RETRIES:
                raise
        sleep(attempt)

    if data is None:
        raise ValueError(f"Unable to load ZFIN data from {url}")

    # Fix weird PMID formatting
    pubs = []
    for publication in data["metaData"]["publications"]:
        pubs.append(publication.replace(": ", ":"))
    data["metaData"]["publications"] = pubs

    # Strip out entries without SO terms:
    valid = []
    for entry in data["data"]:
        if "soTermId" not in entry:
            continue
        valid.append(entry)
    data["data"] = valid
    return data
