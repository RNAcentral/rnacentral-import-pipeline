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
import codecs

import attr
import typing
from pathlib import Path
from pathlib import PurePosixPath

from attr.validators import instance_of as is_a

import textblob as tb

from rnacentral_pipeline.databases.data import IdReference
from . import blob_building


@attr.s()
class WordMatch(object):
    word = attr.ib(validator=is_a(tb.Word))
    name: str = attr.ib(validator=is_a(str))
    group: str = attr.ib(validator=is_a(str))


@attr.s()
class MatchingSentence(object):
    sentence: tb.blob.Sentence = attr.ib(validator=is_a(tb.blob.Sentence))
    matches: typing.List[WordMatch] = attr.ib(validator=is_a(list))
    publication_id = attr.ib(validator=is_a(IdReference))

    def writeables(self, *extra):
        sentence = self.sentence.raw
        sentence = sentence.replace("\n", " ")
        sentence = sentence.replace("  ", " ")
        for match in self.matches:
            data = list(extra)
            rest = [
                self.publication_id.normalized_id,
                match.group,
                match.name,
                match.word,
                sentence,
            ]

            yield data + rest


@attr.s()
class SentenceSelector(object):
    words = attr.ib(type=int, default=200)
    matches = attr.ib(type=int, default=25)

    @classmethod
    def build(cls, word_limit=200, match_limit=25, **kwargs):
        return cls(words=word_limit, matches=match_limit)

    def match(self, matcher, sentence, publication_id):
        if len(sentence.words) > self.words:
            return None
        matches = matcher(sentence)
        if 0 < len(matches) < self.matches:
            return MatchingSentence(sentence, matches, publication_id)


@attr.s()
class PatternMatcher(object):
    group = attr.ib(validator=is_a(str))
    patterns: typing.List[typing.Pattern] = attr.ib()
    pattern: typing.Pattern = attr.ib()

    @classmethod
    def build(cls, group, patterns):
        pattern = re.compile("|".join("(:?%s$)" % p for p in patterns), re.IGNORECASE)

        return cls(
            group=str(group),
            patterns=patterns,
            pattern=pattern,
        )

    def __call__(self, sentence):
        matches = []
        for token in sentence.tokens:
            match = re.match(self.pattern, token)
            if match:
                found = []
                for key, value in match.groupdict().items():
                    if not value:
                        continue
                    found.append(WordMatch(token, key, self.group))
                if not found:
                    raise ValueError("Patterns must be named")
                matches.extend(found)
        return matches


@attr.s()
class NameMatcher(object):
    group = attr.ib(validator=is_a(str))
    names = attr.ib(validator=is_a(set), type=typing.Set[typing.Text])

    @classmethod
    def build(cls, group, names):
        return cls(group=str(group), names=set(names))

    @classmethod
    def from_handle(cls, name, handle):
        return cls.build(name, [n.strip() for n in handle.readlines() if n.strip()])

    def __call__(self, sentence):
        matches = []
        for token in sentence.tokens:
            if token not in self.names:
                continue
            matches.append(WordMatch(token, str(token), self.group))
        return matches


def matches(container, selector, matcher):
    for wrapped in container.blobs():
        for sentence in wrapped.blob.sentences:
            match = selector.match(matcher, sentence, wrapped.pub_id)
            if match:
                yield match


def write_matches(container, selector, matcher, output):
    writer = csv.writer(output)
    for match in matches(container, selector, matcher):
        writer.writerows(match.writeables())


def write_pattern_matches(filename, matcher, output, **kwargs):
    selector = SentenceSelector.build(**kwargs)
    container = blob_building.build(filename)
    write_matches(container, selector, matcher, output)


def write_name_matches(filename, name, name_handle, output, **kwargs):
    matcher = NameMatcher.from_handle(name, name_handle)
    selector = SentenceSelector.build(**kwargs)
    container = blob_building.build(filename)
    write_matches(container, selector, matcher, output)
