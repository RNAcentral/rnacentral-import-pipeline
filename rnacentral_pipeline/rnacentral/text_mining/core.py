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

import six
import attr
from attr.validators import instance_of as is_a
import typing

try:
    from pathlib import Path
    from pathlib import PurePosixPath
except ImportError:
    from pathlib2 import Path
    from pathlib2 import PurePosixPath

import textblob as tb

from . import blob_building


@attr.s()
class WordMatch(object):
    word = attr.ib(type=tb.Word)
    name = attr.ib(type=six.text_type)
    group = attr.ib(validator=is_a(six.text_type))


@attr.s()
class MatchingSentence(object):
    sentence = attr.ib(type=tb.blob.Sentence)
    matches = attr.ib(type=typing.List[WordMatch])
    publication_id = attr.ib(validator=is_a(six.text_type))

    def writeables(self, *extra):
        sentence = self.sentence.raw
        sentence = sentence.replace('\n', ' ')
        sentence = sentence.replace('  ', ' ')
        for match in self.matches:
            data = list(extra)
            data.extend([match.name, match.word, sentence])
            yield data


@attr.s()
class SentenceSelector(object):
    words = attr.ib(type=int, default=200)
    matches = attr.ib(type=int, default=25)

    @classmethod
    def build(cls, word_limit=200, match_limit=25, **kwargs):
        return cls(words=word_limit, matches=match_limit)

    def match(self, matcher, sentence):
        if len(sentence.words) > self.words:
            return None
        matches = matcher(sentence)
        if 0 < len(matches) < self.matches:
            return MatchingSentence(sentence, matches)


@attr.s()
class PatternMatcher(object):
    name = attr.ib(type=six.text_type, validator=is_a(six.text_type))
    patterns = attr.ib(type=typing.List[typing.Pattern])
    pattern = attr.ib(type=typing.Pattern)

    @classmethod
    def build(cls, name, patterns):
        pattern = re.compile(
            '|'.join('(:?%s$)' % p for p in patterns),
            re.IGNORECASE)

        return cls(
            name=name,
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
                    found.append(WordMatch(word=token, name=key))
                if not found:
                    raise ValueError("Patterns must be named")
                matches.extend(found)
        return matches


@attr.s()
class NameMatcher(object):
    names = attr.ib(type=typing.Set[typing.Text])

    @classmethod
    def build(cls, names):
        return cls(names=set(names))

    @classmethod
    def from_handle(cls, handle):
        return cls.build([n.strip() for n in handle.readlines() if n.strip()])

    def __call__(self, sentence):
        return [WordMatch(t, str(t)) for t in sentence.tokens if t in self.names]


def matches(wrapper, selector, matcher):
    for wrapped in container.blobs():
        for sentence in wrapper.blob.sentences:
            match = selector.match(matcher, sentence)
            if match:
                data.append(match)
            yield match


def write_matches(container, selector, matcher, output):
    writer = csv.writer(output)
    basename = PurePosixPath(filename).stem
    for match in matches(container, selector, matcher):
        writer.writerows(match.writeables(blob.pub_id, matcher.name))


def write_matches(filename, selector, matcher, output):
    path = Path(filename)
    container = blob_building.build(path)
    write_matches(container, selector, matcher, output)


def write_pattern_matches(filename, matcher, output, **kwargs):
    selector = SentenceSelector.build(**kwargs)
    write_matches(filename, selector, matcher, output)


def write_name_matches(filename, name_handle, output, **kwargs):
    matcher = NameMatcher.from_handle(name_handle)
    selector = SentenceSelector.build(**kwargs)
    write_matches(filename, selector, matcher, output)
