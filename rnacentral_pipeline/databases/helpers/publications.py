# -*- coding: utf-8 -*-

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

import os
import re
import csv
import json
import logging
import sqlite3
from glob import glob
from xml.etree import cElementTree as ET

import attr
import requests
from retry import retry
from ratelimiter import RateLimiter
from functools32 import lru_cache

from ..data import Reference
from ..data import IdReference


LOGGER = logging.getLogger(__name__)


class UnknownReference(Exception):
    """
    This is created when the requested reference cannot be found either in
    EuropePMC or in the indexed data depending upon how the lookup was done.
    """
    pass


class TooManyPublications(Exception):
    """
    Raised if when we try to lookup a publication in EuropePMC and we find too
    many (> 1) for the given reference.
    """
    pass


TABLE = '''
CREATE TABLE IF NOT EXISTS {name}s (
    id {type} primary key,
    data text,

    UNIQUE (id) ON CONFLICT REPLACE
)
'''


@lru_cache()
@retry(requests.HTTPError, tries=5, delay=1)
@RateLimiter(max_calls=5, period=1)
def summary(id_reference):
    response = requests.get(id_reference.external_url())
    response.raise_for_status()

    data = response.json()
    assert data, "Somehow got no data"
    if data['hitCount'] == 0:
        raise UnknownReference(id_reference)

    if data['hitCount'] > 1:
        raise TooManyPublications(id_reference)

    return data['resultList']['result'][0]


def pretty_location(data):
    issue = data.get('issue', '')
    if issue:
        issue = ('(%s)' % issue)

    pages = data.get('pageInfo', '')
    if 'pageInfo' in data and pages:
        pages = ':' + pages

    location = u'{title} {volume}{issue}{pages} ({year})'.format(
        title=data['journalTitle'],
        issue=issue,
        volume=data.get('journalVolume', ''),
        pages=pages,
        year=data['pubYear'],
    )
    location = location.replace('  ', ' ')
    return location


def reference(ref_id):
    return IdReference.build(ref_id)


def clean_title(title):
    return re.sub(r'\.$', '', title)


def query_pmc(id_reference):
    data = summary(id_reference)
    pmid = data.get('pmid', None)
    if pmid:
        pmid = int(pmid)

    return Reference(
        authors=data['authorString'],
        location=pretty_location(data),
        title=clean_title(data['title']),
        pmid=pmid,
        doi=data.get('doi', None),
    )


def xml_text(tag, node, missing=None, fn=None):
    found = node.find(tag)
    if found is None or not found.text:
        return missing
    if fn:
        return fn(found.text)
    return found.text


def node_to_reference(node):
    pmid = xml_text('pmid', node, fn=int)
    doi = xml_text('DOI', node)
    if not pmid and not doi:
        return None

    authors = []
    for author in node.findall('./AuthorList/Author'):
        last = xml_text('LastName', author)
        initials = xml_text('Initials', author)
        authors.append('%s %s' % (last, initials))
    authors = ', '.join(authors) + '.'

    data = {
        'journalTitle': xml_text('journalTitle', node, missing=''),
        'issue': xml_text('issue', node, missing=''),
        'journalVolume': xml_text('journalVolume', node, missing=''),
        'pageInfo': xml_text('pageInfo', node, missing=''),
        'pubYear': xml_text('pubYear', node, missing=''),
    }

    return Reference(
        authors=authors,
        location=pretty_location(data),
        title=xml_text('title', node, fn=clean_title),
        pmid=pmid,
        doi=doi,
    )


def parse_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    for node in root.findall('./PMC_ARTICLE'):
        ref = node_to_reference(node)
        if not ref:
            continue
        yield ref


def index_xml_directory(directory, output):
    conn = sqlite3.connect(output)
    tables = {'pmid': 'int', 'doi': 'text'}
    for name, pid_type in tables.items():
        conn.execute(TABLE.format(name=name, type=pid_type))

    for filename in glob(os.path.join(directory, '*.xml')):
        with open(filename, 'r') as xml_file:
            for ref in parse_xml(xml_file):
                cursor = conn.cursor()
                for table in tables.keys():
                    key = getattr(ref, table, None)
                    if not key:
                        continue
                    stmt = 'INSERT INTO %ss VALUES(?, ?)' % table
                    data = json.dumps(attr.asdict(ref))
                    cursor.execute(stmt, (key, data))
            conn.commit()
    conn.close()


def query_database(cursor, id_ref, allow_fallback=False):
    query = 'SELECT data from %ss WHERE id=?' % id_ref.namespace
    raw_ref = cursor.execute(query, (id_ref.external_id,)).fetchone()
    if not raw_ref:
        if allow_fallback:
            return query_pmc(id_ref)
        raise UnknownReference(id_ref)
    return Reference(**json.loads(raw_ref[0]))


def write_lookup(db, handle, output, allow_fallback=False):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    reader = csv.reader(handle)
    writer = csv.writer(output)
    for (ref_id, accession) in reader:
        id_ref = IdReference.build(ref_id)
        try:
            ref = query_database(cursor, id_ref, allow_fallback=allow_fallback)
        except UnknownReference:
            LOGGER.warning("Could not find reference for %s", id_ref)
            continue
        writer.writerows(ref.writeable(accession))
    conn.close()


def from_file(handle, output):
    reader = csv.reader(handle)
    writer = csv.writer(output)
    for (ref_id, accession) in reader:
        id_ref = IdReference.build(ref_id)
        try:
            complete = query_pmc(id_ref)
        except UnknownReference:
            LOGGER.warning("Could not find reference for %s", id_ref)
            continue
        writer.writerows(complete.writeable(accession))
