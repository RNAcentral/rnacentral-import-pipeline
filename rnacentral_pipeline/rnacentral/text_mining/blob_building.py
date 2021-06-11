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
import codecs
import itertools as it
from xml.etree import cElementTree as ET

import attr
from attr.validators import instance_of as is_a

from pathlib import Path

import textblob as tb

from rnacentral_pipeline.databases.data import IdReference
from rnacentral_pipeline.databases.helpers import publications as pubs

# original implementation of .itertext() for Python 2.7
# Taken from: https://github.com/python/cpython/blob/2.7/Lib/xml/etree/ElementTree.py#L498
def itertext(node):
    tag = node.tag
    if not isinstance(tag, basestring) and tag is not None:
        return
    if node.text:
        yield node.text
    for e in node:
        for s in itertext(e):
            yield s
        if e.tail:
            yield e.tail


def extract_node_text(node):
    lines = list(node.itertext())
    if not lines:
        return None

    paragraph = lines[0]
    for line in lines[1:]:
        if re.match(r"^[,.)]", line):
            paragraph += line
        else:
            paragraph += " " + line

    if not paragraph:
        return None
    paragraph = paragraph.replace("\n", " ")
    paragraph = re.sub(r"\s+", " ", paragraph)
    return paragraph


def article_as_blob(article):
    text = []
    nodes = it.chain(
        article.findall(".//article-title"),
        article.findall(".//abstract//p"),
        article.findall(".//body//p"),
    )

    for node in nodes:
        paragraph = extract_node_text(node)
        if not paragraph:
            continue
        text.append(paragraph)
    return tb.TextBlob("\n".join(text))


def article_id_reference(article):
    ids = {}
    for pub_id in article.findall("./front/article-meta/article-id"):
        namespace = pub_id.attrib["pub-id-type"]
        ids[namespace] = pub_id.text
    for name in ["pmid", "doi", "pmcid"]:
        if name in ids:
            return pubs.reference("%s:%s" % (name, ids[name]))
    raise ValueError("Could not find an pub id for %s" % ET.tostring(article))


@attr.s()
class TextBlobWrapper(object):
    pub_id = attr.ib(validator=is_a(IdReference))
    blob = attr.ib(validator=is_a(tb.TextBlob))


@attr.s()
class TextBlobContainer(object):
    path = attr.ib(validator=is_a(Path))

    def blobs(self):
        if self.path.is_dir():
            for subpath in self.path.iterdir():
                for blob in self.__blob__(subpath):
                    yield blob

        else:
            for blob in self.__blob__(self.path):
                yield blob

    def __blob__(self, path):
        if path.suffix == ".txt":
            yield self.__text_blob__(path)

        elif path.suffix == ".xml":
            with path.open("r") as raw:
                header = raw.read(10)

            if header.startswith("<articles>"):
                for blob in self.__full_text_blob__(path):
                    yield blob

            elif header.startswith("<PMCSet>"):
                for blob in self.__metadata_blob__(path):
                    yield blob

        else:
            raise ValueError("Cannot parse publication data from: %s" % path)

    def __text_blob__(self, path):
        with codecs.open(path, "r", errors="ignore") as raw:
            text = raw.read()
            blob = tb.TextBlob(text)
        return TextBlobWrapper(pubs.reference(path.stem), blob)

    def __metadata_blob__(cls, path):
        with open(str(path), "r") as raw:
            for ref in pubs.parse_xml(raw):
                id_ref = ref.id_reference()
                blob = tb.TextBlob(ref.title)
                yield TextBlobWrapper(id_ref, blob)

    def __full_text_blob__(self, path):
        with path.open("r") as raw:
            tree = ET.parse(raw)
            root = tree.getroot()
            for node in root.findall("./article"):
                id_ref = article_id_reference(node)
                blob = article_as_blob(node)
                yield TextBlobWrapper(id_ref, blob)


def build(path):
    if isinstance(path, str):
        return TextBlobContainer(Path(path))
    return TextBlobContainer(path)
