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

import itertools as it

import attr

from databases.ena.parsers import parse as ena


def as_entry(entry):
    """
    Modify an ENA entry into an approbate RefSeq entry.
    """
    return attr.evolve(
        entry,
        database='REFSEQ',
        exons=[],
    )


def parse(handle):
    """
    Parse all entries in the handle to produce an iterable of all RefSeq
    entries.
    """
    return it.imap(as_entry, ena(handle))
