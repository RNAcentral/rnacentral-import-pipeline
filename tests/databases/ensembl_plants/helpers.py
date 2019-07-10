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

from rnacentral_pipeline.databases.ensembl_plants import parser


def parse(filename):
    with open(filename, 'r') as raw:
        return list(parser.parse(raw))


def entries_for(entries, accession):
    return [e for e in entries if e.accession == accession]


def entry_for(entries, accession):
    val = entries_for(entries, accession)
    assert len(val) == 1
    return val[0]


def has_entry_for(entries, accession):
    return bool(entries_for(entries, accession))