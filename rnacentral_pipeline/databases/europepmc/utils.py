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

import re
import csv
import logging

from six.moves import html_parser

from rnacentral_pipeline.databases.helpers.publications import reference

LOGGER = logging.getLogger(__name__)


def pretty_location(data):
    """
    Produce the noramlized location description we use from the EuropePMC data.
    """

    issue = data.get("issue", "")
    if issue:
        issue = "(%s)" % issue

    pages = data.get("pageInfo", "")
    if "pageInfo" in data and pages:
        pages = ":" + pages

    location = u"{title} {volume}{issue}{pages} ({year})".format(
        title=data.get("journalTitle", ""),
        issue=issue,
        volume=data.get("journalVolume", ""),
        pages=pages,
        year=data["pubYear"],
    )
    location = location.replace("  ", " ")
    if location.endswith("."):
        return location[0:-1]
    return location


def clean_title(title):
    """
    Cleanup the title into a normalized setup.
    """
    stripped = re.sub(r"\.$", "", title)
    parser = html_parser.HTMLParser()
    return parser.unescape(stripped)


def id_refs_from_handle(handle, column):
    for row in csv.reader(handle):
        raw = row[column]
        rest = [d for i, d in enumerate(row) if i != column]
        yield (reference(raw), rest)


def write_lookup(handler, handle, output, column=0, ignore_errors=False):
    writer = csv.writer(output)
    for id_ref, rest in id_refs_from_handle(handle, column):
        try:
            ref = handler(id_ref)
        except Exception as err:
            LOGGER.warning("Failed to lookup: %s", id_ref)
            LOGGER.exception(err)
            if not ignore_errors:
                raise err
            continue
        writer.writerows(ref.writeable(rest))
