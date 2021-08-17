# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.data.databases import Database

from . import data

REP_DBS = {
    Database.pdbe,
    Database.refseq,
    Database.ensembl,
    Database.hgnc,
}


def should_reject(location: data.LocationInfo) -> bool:
    if any(d in REP_DBS for d in location.databases):
        return False
    return location.qa.has_issue


def should_highlight(location: data.LocationInfo) -> bool:
    qa = location.qa
    if any(d in REP_DBS for d in location.databases):
        return True
    if not qa.has_issue and not location.has_introns():
        return True
    return False


def classify_cluster(state: data.State, cluster: int):
    for location in state.members_of(cluster):
        if should_highlight(location):
            state.highlight_location(location.id)
