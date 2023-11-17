# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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
import json
import operator as op
import typing as ty

from rnacentral_pipeline.databases import data

from . import helpers


def as_expression(mapping):

    pass


def parse(handle, db_url):
    """
    Process the jsonlines output from Rust into entries.

    The jsonlines output is already grouped by geneID and urs taxid so this
    should give us the transcript level linkage we're after without any further
    processing.
    """
    for line in handle:
        hit = json.loads(line)
        for experiment in hit["experiment"]:
            yield helpers.as_entry(hit, experiment)
