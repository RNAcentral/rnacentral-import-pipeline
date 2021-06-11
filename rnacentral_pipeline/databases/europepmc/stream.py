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
import collections as coll

from rnacentral_pipeline.databases.helpers.publications import reference
from rnacentral_pipeline.databases.europepmc import fetch

from . import xml

LOGGER = logging.getLogger(__name__)


def load_ids(handle, column):
    data = coll.defaultdict(list)
    reader = csv.reader(handle)
    for row in reader:
        key = reference(row[column])
        rest = [d for i, d in enumerate(row) if i != column]
        data[key].append(rest)
    return dict(data)


def fallback(data):
    for id_ref, rows in data.items():
        try:
            ref = fetch.lookup(id_ref)
            yield id_ref, ref, rows
        except Exception:
            pass


def lookup(ids, directory, column, allow_fallback=True, ignore_missing=False):
    data = load_ids(ids, column)
    not_found = dict(data)
    for ref in xml.parse_directory(directory):
        for possible_id in ref.id_references:
            if possible_id in data:
                yield ref, data[possible_id]
                del not_found[possible_id]
                if not not_found:
                    return

    if allow_fallback:
        for id_ref, ref, rows in fallback(dict(not_found)):
            yield ref, rows
            del not_found[id_ref]

    for id_ref, rows in not_found.items():
        LOGGER.warning("Failed to lookup: %s (%i entries)", id_ref, len(rows))
        if not ignore_missing:
            raise ValueError("Could not lookup %s" % id_ref)


def write_lookup(
    ids, directory, output, column=0, allow_fallback=True, ignore_missing=False
):
    writer = csv.writer(output)
    for ref, rows in lookup(
        ids,
        directory,
        column,
        allow_fallback=allow_fallback,
        ignore_missing=ignore_missing,
    ):
        for rest in rows:
            writer.writerows(ref.writeable(rest))
