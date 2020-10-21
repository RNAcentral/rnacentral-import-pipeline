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

import csv
import itertools as it
import json
import logging
import operator as op
import typing as ty

from intervaltree import IntervalTree
from more_itertools import partition

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree

from . import data, rrna

LOGGER = logging.getLogger(__name__)

Key = str


def load(handle) -> ty.Iterable[data.UnboundLocation]:
    ontology = so_tree.load_ontology(so_tree.REMOTE_ONTOLOGY)
    for entry in psql.json_handler(handle):
        yield data.UnboundLocation.build(entry, ontology)


def cluster_key(value: data.UnboundLocation) -> Key:
    return value.extent.chromosome


def always_bad_location(location: data.UnboundLocation) -> bool:
    if location.qa.has_issue:
        if location.extent.chromosome == "MT":
            state = location.qa.as_tuple()
            # If the only issue is contamination, it is actually good as we
            # misclassified mito sequences sometimes.
            LOGGER.debug("State: %s", state)
            return state != (True, False, True, False)
        return True
    return False


def always_ignorable_location(location: data.UnboundLocation) -> bool:
    return location.rna_type.insdc in {"antisense_RNA", "lncRNA"}


def cluster(handle) -> ty.Iterable[data.State]:
    entries = load(handle)
    for (chromosome, locations) in it.groupby(entries, cluster_key):
        state = data.State(chromosome=chromosome)
        for location in locations:
            LOGGER.debug("Testing %s", location.name)
            if always_ignorable_location(location):
                LOGGER.debug("Always Ignoring: %s", location.name)
                state.ignore(location)
                continue

            if always_bad_location(location):
                LOGGER.debug("Always Rejected: %s", location.name)
                state.reject(location)
                continue

            locus = data.Locus.singleton(location)
            overlaps = state.overlaps(location)
            if not overlaps:
                LOGGER.debug("Add singleton %s", location.name)
                update = data.StateUpdate(
                    reject=[], ignore=[], new_clusters=[locus], old_clusers=[]
                )
                state.update(update)
            else:
                logger.debug("Merging into %s", [o.data.name for o in overlaps])
                other = []
                for interval in overlaps:
                    other.append(interval.data)
                locus = locus.merge(other)
                update = data.StateUpdate(new_clusters=[locus], old_clusters=[other])
                state.update(update)

        yield state


def split_cluster(cluster: data.Cluster) -> ty.List[data.Cluster]:
    """
    The idea is to split the locus into different locus depending upon RNA type
    and other factors.
    """
    return [cluster]


def handle_rfam_only(cluster) -> data.StateUpdate:
    rfam_only, rest = partition(lambda m: m.info.is_rfam_only(), cluster.members)
    rfam_only = list(rfam_only)
    rest = list(rest)
    if not rfam_only:
        return StateUpdate.no_op()
    ignored = []
    rejected = []
    allowed = [m.info for m in rest]
    for rfam in rfam_only:
        for other in rest:
            if rfam.overlaps(other):
                rejected.append(rfam)
                break
        else:
            allowed.append(rfam.info)
    updated_cluster = Cluster.from_locations(allowed)
    return StateUpdate(
        reject=rejected,
        ignore=ignored,
        new_clusters=[updated_cluster],
        old_clusters=[cluster],
    )


def reject_pseudogenes(cluster, pseudo) -> data.StateUpdate:
    return data.StateUpdate.no_op()


def build(
    clusters: ty.Iterable[data.State], pseudogenes: IntervalTree
) -> ty.Iterable[data.FinalizedState]:
    for cluster in clusters:
        if not cluster.has_clusters():
            yield cluster.finalize()
            continue

        state = data.State(chromosome=cluster.chromosome)
        for cluster in cluster.clusters():
            update = reject_pseudogenes(cluster, pseudo)
            state.update(update)

            for cluster in split_cluster(cluster):
                update = handle_rfam_only(cluster)
                update.validate()
                state.update(update)

                if cluster.rna_type.is_a("rRNA"):
                    update = rrna.classify_cluster(cluster)
                    state.apply_update(update)
                else:
                    state.ignore_all(cluster.locations())
        yield state.finalize()


def from_json(handle) -> ty.Iterable[data.FinalizedState]:
    clusters = cluster(handle)
    return build(clusters, IntervalTree())
