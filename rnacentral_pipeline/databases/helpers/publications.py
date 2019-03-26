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
from glob import glob
from xml.etree import cElementTree as ET

import six

try:
    import dbm
except:
    from six.moves import dbm_gnu as dbm

import attr
from attr.validators import instance_of as is_a

import typing
import requests
from retry import retry
from ratelimiter import RateLimiter

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

from contextlib import contextmanager
try:
    from contextlib import ExitStack
except:
    from contextlib2 import ExitStack

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


@attr.s()
class CacheStorage(object):
    path = attr.ib(validator=is_a(Path))

    def store(self, key, json_data):
        if not key:
            return None
        k = self.__normalize_key__(key)
        self.db[k] = json_data

    def get(self, id_ref, allow_fallback=False):
        key = self.__normalize_key__(id_ref.external_id)
        data = self.db.get(key, None)
        if data:
            return Reference(**json.loads(data))
        if allow_fallback:
            return query_pmc(id_ref)
        raise UnknownReference("Never indexed %s", id_ref)

    @contextmanager
    def open(self, mode='r'):
        self.db = dbm.open(str(self.path), mode)
        yield self
        self.db.close()

    def __normalize_key__(self, raw):
        if isinstance(raw, six.string_types):
            return str(raw.encode('ascii', 'ignore'))
        return str(raw)

@attr.s()
class Cache(object):
    filename = attr.ib(type=Path, validator=is_a(Path))
    keys = attr.ib(
        type=typing.List[six.text_type],
        default=['pmid', 'doi', 'pmcid'],
    )

    @classmethod
    def build(cls, base_path):
        path = Path(base_path)
        if not path.exists():
            path.mkdir()
        return cls(path)

    @classmethod
    def populate_with(cls, base_path, references):
        cache = cls.build(base_path)
        with cache.open(mode='c') as handles:
            for reference in references:
                data = json.dumps(attr.asdict(reference))
                for key_name in cache.keys:
                    key = getattr(reference, key_name, None)
                    handles[key_name].store(key, data)

    @contextmanager
    def open(self, mode='r'):
        with ExitStack() as stack:
            handles = {}
            for key in self.keys:
                filename = self.filename / Path(key)
                cache = CacheStorage(filename)
                db = stack.enter_context(cache.open(mode))
                handles[key] = db
            yield handles


@lru_cache()
@retry(requests.HTTPError, tries=5, delay=1)
@RateLimiter(max_calls=5, period=1)
def summary(id_reference):
    LOGGER.info("Fetching remote summary for %s", id_reference)
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
        authors=six.text_type(data['authorString']),
        location=six.text_type(pretty_location(data)),
        title=six.text_type(clean_title(data['title'])),
        pmid=pmid,
        doi=six.text_type(data.get('doi', None)),
        pmcid=data.get('pmcid', None),
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
    pmcid = six.text_type(xml_text('pmcid', node))

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
        authors=six.text_type(authors),
        location=six.text_type(pretty_location(data)),
        title=six.text_type(xml_text('title', node, fn=clean_title)),
        pmid=pmid,
        doi=six.text_type(doi),
        pmcid=pmcid,
    )


def parse_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    for node in root.findall('./PMC_ARTICLE'):
        ref = node_to_reference(node)
        if not ref:
            continue
        yield ref


def parse_xml_directory(directory):
    for filename in glob(os.path.join(directory, '*.xml')):
        with open(filename, 'r') as xml_file:
            for ref in list(parse_xml(xml_file)):
                yield ref


def index_xml_directory(directory, output):
    Cache.populate_with(output, parse_xml_directory(directory))


def id_refs_from_handle(handle, column=0):
    for row in csv.reader(handle):
        raw = row[column]
        rest = [d for i, d in enumerate(row) if i != column]
        yield (reference(raw), rest)


def write_file_lookup(cache_path, handle, output, column=0, allow_fallback=False):
    writer = csv.writer(output)
    with Cache.build(cache_path).open() as db:
        for id_ref, rest in id_refs_from_handle(handle, column=column):
            try:
                store = db[id_ref.namespace]
                ref = store.get(id_ref, allow_fallback=allow_fallback)
            except UnknownReference:
                LOGGER.warning("Could not handle find reference for %s", id_ref)
                continue
            writer.writerows(ref.writeable(rest))
