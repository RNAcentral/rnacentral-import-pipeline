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

from textblob import TextBlob


def matches(patterns, text):
    pattern = '|'.join('(:?%s$)' % p for p in patterns)
    blob = TextBlob(text.read())
    for word in blob.words:
        if re.match(pattern, word, re.IGNORECASE):
            yield word


def write_matches(patterns, text, output):
    csv.writer(output).writerows(matches(patterns, text))
