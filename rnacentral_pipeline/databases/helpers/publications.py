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

import re
import csv
import json
import sqlite3
from xml.etree import cElementTree as ET

import attr
import requests
from retry import retry
from functools32 import lru_cache

from ..data import Reference
from ..data import IdReference


class UnknownReference(Exception):
    pass


class TooManyPublications(Exception):
    pass


@lru_cache()
@retry(requests.HTTPError, tries=5, delay=1)
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
    print(data)
    issue = data.get('issue', '')
    if issue:
        issue = ('(%s)' % issue)

    pages = data.get('pageInfo', '')
    if 'pageInfo' in data and pages:
        pages = ':' + pages

    location = '{title} {volume}{issue}{pages} ({year})'.format(
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


def lookup_reference(id_reference):
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
        last = author.find('LastName').text
        initials = author.find('Initials').text
        authors.append('%s %s' % (last, initials))
    authors = ', '.join(authors) + '.'

    data = {
        'journalTitle': xml_text('journalTitle', node, missing=''),
        'issue': xml_text('issue', node, missing=''),
        'journalVolumne': xml_text('journalVolume', node, missing=''),
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


def parse_xmls(xml_files):
    for xml_file in xml_files:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for node in root.findall('./PMC_ARTICLE'):
            ref = node_to_reference(node)
            if not ref:
                continue
            yield ref


def index_xml(output, xml_files):
    conn = sqlite3.connect(output)
    tables = {'pmid': 'int', 'doi': 'text'}
    for name, pid_type in tables.items():
        stmt = 'CREATE TABLE {name}s (id {type} primary key, data text)'
        conn.execute(stmt.format(name=name, type=pid_type))

    for ref in parse_xmls(xml_files):
        cursor = conn.cursor()
        for table in tables.keys():
            key = getattr(ref, table, None)
            if not key:
                continue
            stmt = 'INSERT INTO %ss VALUES(?, ?)' % table
            data = json.dumps(attr.asdict(ref))
            cursor.execute(stmt, key, data)


def lookup_refs(db, handle):
    conn = sqlite3.connect(db)
    reader = csv.reader(handle)
    for (ref_id, accession) in reader:
        cursor = conn.cursor()
        id_ref = IdReference.build(ref_id)
        query = 'SELECT data from %ss WHERE id=?' % id_ref.namespace
        raw_ref = cursor.execute(query, id_ref.external_id).fetchone()
        if not raw_ref:
            raise UnknownReference(id_ref)
        ref = Reference(**json.loads(raw_ref))
        yield accession, ref


def write_lookup(db, handle, output):
    writer = csv.writer(output)
    for accession, ref in lookup_refs(db, handle):
        writer.writerows(ref.writeable(accession))


def from_file(handle, output):
    reader = csv.reader(handle)
    writer = csv.writer(output)
    for (ref_id, accession) in reader:
        id_ref = IdReference.build(ref_id)
        complete = lookup_reference(id_ref)
        writer.writerows(complete.writeable(accession))
