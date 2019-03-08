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

import six
import attr
import typing

import textblob as tb


@attr.s()
class MatchingSentence(object):
    sentence = attr.ib(type=tb.blob.Sentence)
    matches = attr.ib(type=typing.List[tb.Word])

    def writeables(self, *extra):
        for match in self.matches:
            data = list(extra)
            data.extend([self.filename, match])
            yield data


def select_sentences(pattern, blob, word_limit=200, match_limit=25):
    for sentence in blob.sentences:
        if len(sentence.words) > word_limit:
            continue
        matches = [t for t in sentence.tokens if re.match(pattern, t)]
        if 0 < len(matches) < match_limit:
            yield MatchingSentence(sentence, matches)


def matches(patterns, text):
    pattern = '|'.join('(:?%s$)' % p for p in patterns)
    blob = tb.TextBlob(text)
    return select_sentences(pattern, blob)
