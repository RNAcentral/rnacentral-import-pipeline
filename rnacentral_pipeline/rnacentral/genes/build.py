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
from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree
from rnacentral_pipeline.rnacentral.genes import data, rrna

LOGGER = logging.getLogger(__name__)


def load(handle) -> ty.Iterable[data.LocationInfo]:
    ontology = so_tree.load_ontology(so_tree.REMOTE_ONTOLOGY)
    for entry in psql.json_handler(handle):
        yield data.LocationInfo.build(entry, ontology)


def always_bad_location(location: data.LocationInfo) -> bool:
    if location.qa.has_issue:
        if location.extent.chromosome == "MT":
            state = location.qa.as_tuple()
            # If the only issue is contamination, it is actually good as we
            # misclassified mito sequences sometimes.
            LOGGER.debug("State: %s", state)
            return state != (True, False, True, False)
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
    for cluster in clusters:
        if not cluster.as_interval().overlaps(location.as_interval()):
            continue
        if not has_compatible_rna_types(location, cluster):
            continue
        to_merge.append(cluster)

    if not to_merge:
        return None
    return to_merge


def handle_rfam_only(state: data.State, cluster: int):
    members = state.members_of(cluster)
    if len(members) == 1:
        return

    for member in members:
        if member.location.is_rfam_only():
            state.reject_location(member.location)


def overlaps_pseudogene(location: data.LocationInfo, pseudo: IntervalTree) -> bool:
    return pseudo.overlaps(location.as_interval())


def build(
    locations: ty.Iterable[data.LocationInfo], pseudogenes: IntervalTree
) -> ty.Iterable[data.FinalizedState]:
    for (key, locations) in it.groupby(locations, data.ClusteringKey.from_location):
        state = data.State(key=key)
        for location in locations:
            state.add_location(location)
            LOGGER.debug("Testing %s", location.name)

            if always_ignorable_location(location):
                LOGGER.debug("Always Ignoring: %s", location.name)
                state.ignore_location(location)
                continue

            if always_bad_location(location):
                LOGGER.debug("Always Rejected: %s", location.name)
                state.reject_location(location)
                continue

            if overlaps_pseudogene(location, pseudogenes):
                LOGGER.debug("Rejecting Pseudogene: %s", location.name)
                state.reject_location(location)
                continue

            overlaps = state.overlaps(location)
            if not overlaps:
                LOGGER.debug("Adding singleton cluster of %s", location)
                state.add_singleton_cluster(location)
                continue

            possible = [i.data for i in overlaps]
            to_merge = select_mergable(location, possible)
            if to_merge is None:
                LOGGER.debug("Adding singleton cluster of %s", location)
                state.add_singleton_cluster(location)
                continue

            LOGGER.debug("Merging into %s", [c.name for c in to_merge])
            state.merge_clusters(to_merge, additional=[location])

        if not state.has_clusters():
            yield state.finalize()
            continue

        for cluster_id in state.clusters():
            handle_rfam_only(state, cluster_id)

            members = state.members_of(cluster_id)
            if any(m.location.rna_type.is_a("rRNA") for m in members):
                rrna.classify_cluster(state, cluster_id)
            else:
                state.ignore_cluster(cluster_id)
        yield state.finalize()


def from_json(handle) -> ty.Iterable[data.FinalizedState]:
    locations = load(handle)
    return build(locations, IntervalTree())
