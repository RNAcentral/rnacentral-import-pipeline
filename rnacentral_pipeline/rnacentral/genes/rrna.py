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

from . import data

REP_DBS = {
    "pdbe",
    "refseq",
    "ensembl",
    "hgnc",
}


def should_reject(location: data.LocationInfo) -> bool:
    if any(d.lower() in REP_DBS for d in location.databases):
        return False
    return location.qa.has_issue


def should_highlight(member: data.ClusterMember) -> bool:
    qa = member.location.qa
    if any(d.lower() in REP_DBS for d in member.location.databases):
        return True
    if not qa.has_issue and not member.location.has_introns():
        return True
    return False


def classify_cluster(state: data.State, cluster: int):
    for member in state.members_of(cluster):
        if should_highlight(member):
            state.highlight_location(member.location.id)
