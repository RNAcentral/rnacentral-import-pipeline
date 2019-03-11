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

from . import core


def matches(patterns, text, word_limit=200, match_limit=25):
    pattern = '|'.join('(:?%s$)' % p for p in patterns)
    blob = tb.TextBlob(text)
    def fn(sent):
        return [t for t in sent.tokens if re.match(pattern, t, re.IGNORECASE)]

    selector = core.SentenceSelector(words=word_limit, matches=match_limit)
    for match in core.select_sentences(blob, selector, fn):
        yield match
