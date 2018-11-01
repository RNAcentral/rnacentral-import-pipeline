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
from time import sleep
from xml.etree import cElementTree as ET

import attr
import requests
from functools32 import lru_cache

from ..data import Reference

URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}+AND+SRC:MED&format=json'


class UnknownPmid(Exception):
    pass


class FailedPublicationId(Exception):
    pass


class TooManyPublications(Exception):
    pass


@lru_cache()
def summary(pmid):
    for count in xrange(5):
        response = requests.get(URL.format(pmid=pmid))
        try:
            response.raise_for_status()
            break
        except requests.HTTPError:
            if response.status_code == 500:
                sleep(0.1 * (count + 1))
                continue
        else:
            raise FailedPublicationId(pmid)

    data = response.json()
    assert data, "Somehow got no data"
    if data['hitCount'] == 0:
        raise UnknownPmid(pmid)

    if data['hitCount'] == 1:
        return data['resultList']['result'][0]

    raise TooManyPublications(pmid)


def pretty_location(data):
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


def reference(pmid):
    data = summary(pmid)
    title = re.sub(r'\.$', '', data['title'])
    return Reference(
        authors=data['authorString'],
        location=pretty_location(data),
        title=title,
        pmid=int(pmid),
        doi=data.get('doi', None),
    )


def node_to_reference(node):
    pmid = node.find('./pmid')
    doi = node.find('./DOI')
    if not pmid and not doi:
        return None

    if pmid:
        pmid = int(pmid.text)

    if doi:
        doi = doi.text

    data = {}
    return Reference(
        authors=authors,
        location=pretty_location(data),
        title=node.find('./title').text,
        pmid=pmid,
        doi=doi,
    )


def index_xml(output, xml_files):
    conn = sqlite3.connect(output)
    conn.execute('CREATE TABLE publications (id int, data text)')
    for xml_file in xml_files:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for node in root.findall('./PMC_ARTICLE'):
            data = node_to_reference(node)
            if not data:
                continue
            pmid = data['pmid']
            data = json.dumps(attr.asdict(data))
            conn.execute('INSERT INTO publications VALUES(?, ?)', pmid, data)


def lookup_refs(db, handle):
    conn = sqlite3.connect(db)
    reader = csv.reader(handle)
    for (pmid, accession) in reader:
        raw_ref = conn.execute(
            "SELECT data from publications WHERE id=?",
            pmid
        )
        ref = Reference(**json.loads(raw_ref))
        yield accession, ref


def write_lookup(db, handle, output):
    writer = csv.writer(output)
    for accession, ref in lookup_refs(db, handle):
        writer.writerows(ref.writeable(accession))


def from_file(handle, output):
    reader = csv.reader(handle)
    writer = csv.writer(output)
    for (pmid, accession) in reader:
        ref = reference(int(pmid))
        writer.writerows(ref.writeable(accession))
