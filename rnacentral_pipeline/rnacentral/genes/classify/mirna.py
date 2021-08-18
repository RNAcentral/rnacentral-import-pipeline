# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import typing as ty

from rnacentral_pipeline.databases.data.databases import Database
from rnacentral_pipeline.rnacentral.genes import data

DBS = [
    Database.mirbase,
    Database.mirgenedb,
]


def ordered_choice_highlight(state: data.State, cluster: int, selector, choices: ty.List[Database]):
    grouped = [[] for x in range(len(choices))]
    for location in state.members_of(cluster):
        if not selector(location):
            continue
        for index, db in enumerate(choices):
            if db in location.databases:
                grouped[index].append(location)

    for locations in grouped:
        if not locations:
            continue
        for location in locations:
            state.highlight_location(location.id)
        return


def is_mature(location: data.LocationInfo) -> bool:
    return location.rna_type.is_a('SO:0001244')


def is_precursor(location: data.LocationInfo) -> bool:
    return location.rna_type.is_a('SO:0000276')


def classify_cluster(state: data.State, cluster: int):
    ordered_choice_highlight(state, cluster, is_precursor, DBS)
    ordered_choice_highlight(state, cluster, is_mature, DBS)
