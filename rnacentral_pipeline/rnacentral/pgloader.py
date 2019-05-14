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


SUCCESS_PATTERN = re.compile(r'^\s+Total import time\s+[?✓]')
FAILED_PATTERN = re.compile(r'^\s+Total import time\s+\d+')


def summary_line(handle):
    summary = None
    for line in handle:
        if line.strip().startswith('Total import time'):
            if summary:
                raise ValueError("Mutliple summaries found")
            summary = line
    if not summary:
        raise ValueError("Did not find a summary line to analyze")
    return summary


def validate(output):
    summary = summary_line(output)
    if re.match(SUCCESS_PATTERN, summary):
        return True
    if re.match(FAILED_PATTERN, summary):
        return False
    raise ValueError("Could not parse summary line: '%s'" % summary)
