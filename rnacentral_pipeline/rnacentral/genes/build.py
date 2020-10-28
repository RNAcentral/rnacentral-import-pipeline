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
        location = data.LocationInfo.build(entry, ontology)
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


def overlaps_pseudogene(location: data.LocationInfo, pseudo: IntervalTree) -> bool:
    LOGGER.debug("Checking %s for overlaps to pseudogenes", location.id)
    overlaps = pseudo.overlaps(location.as_interval())
    if overlaps:
        LOGGER.debug("Overlaps %s", overlaps)
        return True
    LOGGER.debug("None found")
    return False


def build(
    locations: ty.Iterable[data.LocationInfo], pseudogenes: IntervalTree
) -> ty.Iterable[data.FinalizedState]:
    for (key, locations) in it.groupby(locations, data.ClusteringKey.from_location):
        LOGGER.debug("Building clusters for %s", key)
        state = data.State(key=key)
        for location in locations:
            LOGGER.debug("Testing %s", location.id)
            state.add_location(location)
            LOGGER.debug(
                "Lengths, locations: %i, tree: %i, clusters: %i", *state.lengths()
            )

            if always_ignorable_location(location):
                LOGGER.debug("Always Ignoring: %s", location.id)
                state.ignore_location(location)
                continue

            if always_bad_location(location):
                LOGGER.debug("Always Rejected: %s", location.id)
                state.reject_location(location)
                continue

            if overlaps_pseudogene(location, pseudogenes):
                LOGGER.debug("Rejecting Pseudogene: %s", location.id)
                state.reject_location(location)
                continue

            overlaps = state.overlaps(location)
            if not overlaps:
                LOGGER.debug("Adding singleton cluster of %s", location.id)
                state.add_singleton_cluster(location)
                continue

            to_merge = select_mergable(location, overlaps)
            if to_merge is None:
                LOGGER.debug("Adding singleton cluster of %s", location.id)
                state.add_singleton_cluster(location)
                continue
            elif len(to_merge) == 1:
                cluster = to_merge[0]
                LOGGER.debug(
                    "Adding location %i to cluster %i", location.id, cluster.id
                )
                state.add_to_cluster(location, cluster.id)
            else:
                LOGGER.debug("Merging into %s", [c.id for c in to_merge])
                cluster_id = state.merge_clusters(to_merge)
                state.add_to_cluster(location, cluster_id)

        state.validate()
        LOGGER.debug("Done building clusters for %s", key)
        if not state.has_clusters():
            LOGGER.debug("No clusters to analyze")
            yield state.finalize()
            continue

        for cluster_id in state.clusters():
            LOGGER.debug("Analyzing cluster %i", cluster_id)
            handle_rfam_only(state, cluster_id)

            if not state.has_cluster(cluster_id):
                continue

            members = state.members_of(cluster_id)
            if any(l.rna_type.is_a("rRNA") for l in members):
                rrna.classify_cluster(state, cluster_id)
            else:
                state.ignore_cluster(cluster_id)
        yield state.finalize()


def from_json(handle) -> ty.Iterable[data.FinalizedState]:
    locations = load(handle)
    return build(locations, IntervalTree())
