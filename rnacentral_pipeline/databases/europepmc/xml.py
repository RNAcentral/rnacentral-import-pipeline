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

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from contextlib import contextmanager
try:
    from contextlib import ExitStack
except:
    from contextlib2 import ExitStack

from .utils import clean_title
from .utils import pretty_location

from ..data import Reference


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
        assert key, "Provided no key for index of: %s" % json_data
        assert json_data
        k = self.__normalize_key__(key)
        assert k, "Failed to get normalized key of %s" % key
        self.db[k] = json_data

    def get(self, id_ref, allow_fallback=False):
        key = self.__normalize_key__(id_ref.external_id)
        data = self.db.get(key, None)
        if data:
            return Reference(**json.loads(data))
        if allow_fallback:
            return query_pmc(id_ref)
        raise UnknownReference("Never indexed %s" % id_ref)

    @contextmanager
    def open(self, mode='r'):
        self.db = dbm.open(str(self.path), mode)
        yield self
        self.db.close()

    def __normalize_key__(self, raw):
        if isinstance(raw, six.string_types):
            return raw.encode('ascii', 'ignore')
        return str(raw)

@attr.s()
class Cache(object):
    filename = attr.ib(type=Path, validator=is_a(Path))

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
                for key_name in KnownServices:
                    key = getattr(reference, key_name, None)
                    if key:
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


def write_file_lookup(cache_path, handle, output, column=0, allow_fallback=False,
                      ignore_missing=True):
    writer = csv.writer(output)
    with Cache.build(cache_path).open() as db:
        for id_ref, rest in id_refs_from_handle(handle, column=column):
            try:
                store = db[id_ref.namespace]
                ref = store.get(id_ref, allow_fallback=allow_fallback)
            except UnknownReference as err:
                LOGGER.warning("Could not handle find reference for %s", id_ref)
                if not ignore_missing:
                    raise err
                continue
            except TooManyPublications as err:
                LOGGER.warning("Found too many publications for %s", id_ref)
                if not ignore_missing:
                    raise err
                continue
            writer.writerows(ref.writeable(rest))
