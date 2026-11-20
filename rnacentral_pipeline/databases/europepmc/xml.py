# -*- coding: utf-8 -*-

from __future__ import absolute_import

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import logging
import os
from glob import glob

from lxml import etree as ET

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.europepmc.utils import (
    clean_title,
    pretty_location,
)

LOGGER = logging.getLogger(__name__)

COMMIT_SIZE = 100000


def xml_text(tag, node, missing=None, fn=None):
    found = node.find(tag)
    if found is None or not found.text:
        return missing
    if fn:
        return fn(found.text)
    return found.text


def node_to_reference(node):
    pmid = xml_text("pmid", node, fn=int)
    doi = xml_text("DOI", node)
    if not pmid and not doi:
        return None
    pmcid = str(xml_text("pmcid", node))

    authors = []
    for author in node.findall("./AuthorList/Author"):
        last = xml_text("LastName", author)
        initials = xml_text("Initials", author)
        authors.append("%s %s" % (last, initials))
    authors = ", ".join(authors) + "."

    data = {
        "journalTitle": xml_text("journalTitle", node, missing=""),
        "issue": xml_text("issue", node, missing=""),
        "journalVolume": xml_text("journalVolume", node, missing=""),
        "pageInfo": xml_text("pageInfo", node, missing=""),
        "pubYear": xml_text("pubYear", node, missing=""),
    }

    return Reference(
        authors=str(authors),
        location=str(pretty_location(data)),
        title=str(xml_text("title", node, fn=clean_title)),
        pmid=pmid,
        doi=str(doi),
        pmcid=pmcid,
    )


def parse(xml_file):
    for _, node in ET.iterparse(
        xml_file, recover=True, events=("end",), tag="PMC_ARTICLE"
    ):
        ref = node_to_reference(node)
        if not ref:
            continue
        yield ref
        # Prevents building the whole XML tree in memory.
        # https://stackoverflow.com/questions/9809469/python-sax-to-lxml-for-80gb-xml/9814580#9814580
        # Could remove more things but this is ok.
        node.clear()


def parse_directory(directory):
    for filename in glob(os.path.join(directory, "*.xml")):
        for ref in parse(filename):
            yield ref
