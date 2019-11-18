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

import os
import re
import csv
import dbm
import json
import typing
import logging
from glob import glob
import functools as ft
from pathlib import Path
from contextlib import ExitStack
from contextlib import contextmanager
from xml.etree import cElementTree as ET

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.data import IdReference
from rnacentral_pipeline.databases.data import KnownServices
from rnacentral_pipeline.databases.helpers.publications import reference

from rnacentral_pipeline.databases.europepmc.utils import clean_title
from rnacentral_pipeline.databases.europepmc.utils import pretty_location
from rnacentral_pipeline.databases.europepmc.utils import write_lookup
from rnacentral_pipeline.databases.europepmc.fetch import lookup


LOGGER = logging.getLogger(__name__)


class UncachedReference(Exception):
    """
    This is created when the requested reference was never indexed. This may be
    because it is not an open access publication, or it does not exist at all.
    """
    pass


@attr.s()
class Storage(object):
    namespace = attr.ib(validator=is_a(KnownServices))
    path = attr.ib(validator=is_a(Path))

    @classmethod
    def build(cls, service, base_path):
        return Storage(service, base_path / service.name)

    def store(self, id_ref, json_data):
        assert id_ref.namespace is self.namespace
        key = self.__normalize_key__(id_ref.external_id)
        self.db[key] = json_data

    def get(self, id_ref):
        assert id_ref.namespace is self.namespace
        key = self.__normalize_key__(id_ref.external_id)
        data = self.db.get(key, None)
        if not data:
            return None
        return Reference(**json.loads(data))

    @contextmanager
    def open(self, mode='r'):
        self.db = dbm.open(str(self.path), mode)
        yield self
        self.db.close()

    def __normalize_key__(self, raw):
        if isinstance(raw, str):
            encoded = raw.encode('ascii', 'ignore')
        else:
            encoded = str(raw)
        assert encoded
        return encoded

@attr.s()
class Cache(object):
    _path = attr.ib(type=Path, validator=is_a(Path))
    _pmid = attr.ib(validator=is_a(Storage))
    _doi = attr.ib(validator=is_a(Storage))
    _pmcid = attr.ib(validator=is_a(Storage))

    @classmethod
    def build(cls, base_path):
        path = Path(base_path)
        if not path.exists():
            path.mkdir()
        kwargs = {}
        for service in KnownServices:
            kwargs[service.name] = Storage.build(service, path)
        return cls(path=path, **kwargs)

    @classmethod
    def populate_with(cls, base_path, references):
        cache = cls.build(base_path)
        with cache.open(mode='c') as cache:
            for reference in references:
                cache.store(reference)

    @contextmanager
    def open(self, mode='r'):
        with ExitStack() as stack:
            for service in KnownServices:
                storage = self.__storage__(service)
                db = stack.enter_context(storage.open(mode))
            yield self

    def get(self, id_reference, allow_fallback=False):
        storage = self.__storage__(id_reference.namespace)
        data = storage.get(id_reference)
        if data:
            return data
        elif allow_fallback:
            return lookup(id_reference)
        raise UncachedReference(id_reference)

    def store(self, reference):
        data = json.dumps(attr.asdict(reference))
        for id_reference in reference.id_references:
            storage = self.__storage__(id_reference)
            storage.store(id_reference, data)

    def __storage__(self, ref):
        if isinstance(ref, KnownServices):
            name = ref.name
        elif isinstance(ref, IdReference):
            name = ref.namespace.name
        else:
            raise ValueError("Unknown type: %s" % ref)

        key = '_' + name
        if not hasattr(self, key):
            raise ValueError("Never cached references for: %s" % name)
        return getattr(self, key)


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
    pmcid = str(xml_text('pmcid', node))

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
        authors=str(authors),
        location=str(pretty_location(data)),
        title=str(xml_text('title', node, fn=clean_title)),
        pmid=pmid,
        doi=str(doi),
        pmcid=pmcid,
    )


def parse(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    for node in root.findall('./PMC_ARTICLE'):
        ref = node_to_reference(node)
        if not ref:
            continue
        yield ref


def parse_directory(directory):
    for filename in glob(os.path.join(directory, '*.xml')):
        with open(filename, 'r') as xml_file:
            for ref in list(parse(xml_file)):
                yield ref


def index_directory(directory, output):
    Cache.populate_with(output, parse_directory(directory))


def id_refs_from_handle(handle, column=0):
    for row in csv.reader(handle):
        raw = row[column]
        rest = [d for i, d in enumerate(row) if i != column]
        yield (reference(raw), rest)


def write_file_lookup(cache_path, handle, output, column=0, allow_fallback=False,
                      ignore_missing=True):
    with Cache.build(cache_path).open() as db:
        write_lookup(
            ft.partial(db.get, allow_fallback=allow_fallback),
            handle,
            output,
            column=column,
            ignore_errors=ignore_missing,
        )
