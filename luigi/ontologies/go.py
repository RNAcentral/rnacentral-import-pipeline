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

import gzip
import itertools as it

import attr
from attr.validators import instance_of as is_a

from databases.rfam import utils as rfam
from databases.quickgo import parser as quickgo

from . import helpers as ont


@attr.s()
class TermSources(object):
    rfam_version = attr.ib(validator=is_a(basestring))
    quickgo_file = attr.ib(validator=is_a(basestring))


def rfam_ids(source):
    for family in rfam.load_families(source.rfam_version):
        for (go_term_id, _) in family.go_terms:
            yield go_term_id


def quickgo_ids(source):
    with gzip.open(source.quickgo_file, 'r') as handle:
        for annotation in quickgo.parser(handle):
            yield annotation.term_id


def to_load(source):
    go_ids = it.chain(rfam_ids(source), quickgo_ids(source))
    go_ids = {go_id for go_id in go_ids}
    for go_id in go_ids:
        yield ont.term(go_id)
