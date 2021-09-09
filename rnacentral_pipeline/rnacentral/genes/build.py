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

import itertools as it
import logging
import typing as ty

from intervaltree import IntervalTree

from rnacentral_pipeline import psql
from rnacentral_pipeline.rnacentral.genes import data
from rnacentral_pipeline.rnacentral.genes import classify
from rnacentral_pipeline.rnacentral.genes.data import Methods

LOGGER = logging.getLogger(__name__)


def load(context: data.Context, handle: ty.IO) -> ty.Iterable[data.LocationInfo]:
    for entry in psql.json_handler(handle):
        location = data.LocationInfo.build(context, entry)
        yield location


def always_bad_location(location: data.LocationInfo) -> bool:
    if location.qa.has_issue:
        if location.extent.chromosome == "MT":
            state = location.qa.as_tuple()
            # If the only issue is contamination, it is actually good as we
            # misclassified mito sequences sometimes.
            LOGGER.debug("State: %s", state)
            return state != (True, False, True, False)
        return True
    if location.rna_type.is_a('piRNA'):
        return True
    return False


def always_ignorable_location(location: data.LocationInfo) -> bool:
    return False


def has_compatible_rna_types(
    location: data.LocationInfo, cluster: data.Cluster
) -> bool:
    rna_types = cluster.rna_types()
    if not rna_types:
        raise ValueError("This should be impossible")

    if len(rna_types) != 1:
        return False
    rna_type = rna_types.pop()
    return rna_type.normalized_term == location.rna_type.normalized_term


def select_mergable(
    location: data.LocationInfo, clusters: ty.List[data.Cluster]
) -> ty.Optional[ty.List[data.Cluster]]:
    to_merge = []
    target = location.as_interval()
    for cluster in clusters:
        if not cluster.as_interval().overlaps(target):
            continue
        if not has_compatible_rna_types(location, cluster):
            continue
        to_merge.append(cluster)

    if not to_merge:
        return None
    return to_merge


def handle_rfam_only(state: data.State, cluster: int):
    LOGGER.debug("Checking cluster %i for Rfam only exclusions", cluster)
    members = state.members_of(cluster)
    if len(members) == 1:
        LOGGER.debug("Singleton cluster has no Rfam")
        return

    for location in members:
        LOGGER.debug("Checking %s for Rfam exclusions", location.id)
        if location.is_rfam_only():
            LOGGER.debug("Rejecting %s as it is Rfam only", location.id)
            state.reject_location(location)


def overlaps_pseudogene(context: data.Context, location: data.LocationInfo) -> bool:
    LOGGER.debug("Checking %s for overlaps to pseudogenes", location.id)
    overlaps = context.overlaps_pseudogene(location)
    if overlaps:
        LOGGER.debug("Overlaps %s", overlaps)
        return True
    LOGGER.debug("None found")
    return False


def is_pirna_singleton(state: data.State, cluster: int) -> bool:
    members = state.members_of(cluster)
    target = members[0]
    return len(members) == 1 and target.rna_type.is_a('piRNA') and not bool(target.providing_databases)


def build(
        context: data.Context, method: Methods, locations: ty.Iterable[data.LocationInfo]
) -> ty.Iterable[data.FinalizedState]:
    handler = method.handler()
    for (key, locations) in it.groupby(locations, data.ClusteringKey.from_location):
        LOGGER.debug("Building clusters for %s", key)
        state = data.State(key=key, method=method.name)
        for location in locations:
            LOGGER.debug("Testing %s", location.id)
            state.add_location(location)
            LOGGER.debug(
                "Lengths, locations: %i, tree: %i, clusters: %i", *state.lengths()
            )

            handler.handle_location(state, context, location)

        state.validate()
        LOGGER.debug("Done building clusters for %s", key)
        if not state.has_clusters():
            LOGGER.debug("No clusters to analyze")
            yield state.finalize()
            continue

        for cluster_id in state.clusters():
            LOGGER.debug("Analyzing cluster %i", cluster_id)
            if not state.has_cluster(cluster_id):
                raise ValueError("Seems link this should be impossible")

            handle_rfam_only(state, cluster_id)

            if is_pirna_singleton(state, cluster_id):
                LOGGER.debug("Rejecting a piRNA singleton %s", cluster_id)
                for location in state.members_of(cluster_id):
                    state.reject_location(location)


            members = state.members_of(cluster_id)
            if any(l.rna_type.is_a("rRNA") for l in members):
                classify.rrna(state, cluster_id)
            elif any(l.rna_type.is_a('miRNA') for l in members):
                classify.mirna(state, cluster_id)
            else:
                LOGGER.info("Do not know how to cluster cluster %s", cluster_id)
        yield state.finalize()


def from_json(context: data.Context, method: Methods, handle: ty.IO) -> ty.Iterable[data.FinalizedState]:
    locations = load(context, handle)
    return build(context, method, locations)
