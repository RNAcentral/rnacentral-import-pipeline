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

import typing as ty

import attr
from attr.validators import instance_of as is_a
from intervaltree import IntervalTree

from . import data

REP_DBS = {
    "pdbe",
    "refseq",
    "ensembl",
    "hgnc",
}


def should_reject(location: data.UnboundLocation) -> bool:
    if any(d.lower() in REP_DBS for d in location.databases):
        return False
    return location.qa.has_issue


def highlight_members(
    members: ty.List[data.ClusterMember],
) -> ty.List[data.ClusterMember]:
    updates = []
    for member in members:
        qa = member.info.qa
        member_type = data.MemberType.member
        if any(d.lower() in REP_DBS for d in member.info.databases):
            member_type = data.MemberType.highlighted
        elif not qa.has_issue and not member.info.has_introns():
            member_type = data.MemberType.highlighted
        updates.append(attr.assoc(member, member_type=member_type))
    return updates


def classify(cluster: data.Cluster) -> ty.List[data.Cluster]:
    cluster = attr.assoc(cluster, members=highlight_members(cluster.members))
    return [cluster]
