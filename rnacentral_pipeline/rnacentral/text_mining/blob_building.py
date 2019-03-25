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
import itertools as it
from xml.etree import cElementTree as ET

import six
import attr
from attr.validators import instance_of as is_a

try:
    from pathlib import Path
    from pathlib import PurePosixPath
except ImportError:
    from pathlib2 import Path
    from pathlib2 import PurePosixPath

import textblob as tb

from rnacentral_pipeline.databases.helpers import publications as pubs


def article_as_blob(article):
    text = []
    nodes = it.chain(
        article.findall('.//abstract//p'),
        article.findall('.//body//p'),
    )

    for node in nodes:
        paragraph = node.text
        if not paragraph:
            continue
        paragraph = paragraph.replace('\n' ' ')
        paragraph = re.sub(r'\s+', ' ', paragraph)
        text.append(paragraph)
    return tb.TextBlob('\n'.join(text))


def article_id_reference(article):
    ids = {}
    for pub_id in article.findall('article-id'):
        ids[pub_id.attrib['pub-id-type']] = pub_id.text
    for name in ['pmid', 'doi', 'pmcid']:
        if name in ids:
            return pubs.reference('%s:%s' % (name, ids[name]))
    raise ValueError("Could not find an pub id for %s" % ids)


@attr.s()
class TextBlobWrapper(object):
    pub_id = attr.ib(validator=is_a(six.text_type))
    blob = attr.ib(validator=is_a(tb.TextBlob))


@attr.s()
class TextBlobContainer(object):
    path = attr.ib(validator=is_a(Path))

    def blobs(self):
        if self.path.is_dir():
            for subpath in path.iterdir():
                for blob in self.__blob__(subpath):
                    yield blob
        else:
            for blob in self.__blob__(self.path):
                yield blob

    def __blob__(self, path):
        if path.endswith('.txt'):
            yield self.__text_blob__(path)

        if path.endswith('.xml'):
            with path.open('r') as raw:
                header = raw.read(10)
                if header.startswith('<articles>'):
                    yield self.__full_text_blob__(path)
                elif header.startswith('<PMCSet>'):
                    yield self.__metadata_blob__(path)

        else:
            raise ValueError("Cannot parse publication data from: %s" % path)

    def __text_blob__(self, path):
        with codecs.open(path, 'r', errors='ignore') as text:
            blob = tb.TextBlob(text.read().decode('ascii', 'ignore'))
            yield TextBlobWrapper(path.stem, blob)

    def __metadata_blob__(cls, path):
        with path.open('r') as raw:
            for ref in pubs.parse_xml(raw):
                id_ref = ref.id_reference()
                blob = tb.TextBlob(reference.title)
                yield TextBlobWrapper(id_ref.normalized_id(), blob)

    def __full_text_blob__(self, path):
        with path.open('r') as raw:
            tree = ET.parse(raw)
            root = tree.getroot()
            for node in root.findall('./article'):
                id_ref = article_id_reference(node)
                blob = article_as_blob(node)
                yield TextBlobWrapper(id_ref.normalized_id(), blob)


def build(path):
    TextBlobWrapper(path)
