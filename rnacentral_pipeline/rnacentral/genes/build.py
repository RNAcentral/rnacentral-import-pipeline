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
import operator as op
import typing as ty

from intervaltree import IntervalTree
from more_itertools import partition

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree

from . import data, rrna

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
            if always_ignorable_location(location):
                state.ignore(location)
                continue

            if always_bad_location(location):
                state.reject(location)
                continue

            locus = data.Locus.singleton(location)
            overlaps = state.overlaps(location)
            if not overlaps:
                state.add(locus)
            else:
                other = []
                for interval in overlaps:
                    state.remove_interval(interval)
                    other.append(interval.data)
                locus = locus.merge(other)
                state.add(locus)

        yield state


def split_cluster(cluster):
    """
    The idea is to split the locus into different locus depending upon RNA type
    and other factors.
    """
    pass


def handle_rfam_only(cluster):
    rfam_only, rest = partition(lambda m: m.info.is_rfam_only(), cluster.members)
    rfam_only = list(rfam_only)
    rest = list(rest)
    if not rfam_only:
        return (cluster, [], [])
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
    return (updated_cluster, ignored, rejected)


def build(
    clusters: ty.Iterable[data.State], pseudogenes: IntervalTree
) -> ty.Iterable[data.FinalizedState]:
    for cluster in clusters:
        state = data.State(chromosome=cluster.chromosome)
        for cluster in cluster.clusters():
            cluster_interval = cluster.as_interval()
            overlaping_pseudo = psuedogenes.overlaps(cluster_interval)
            if overlaping_pseudo:
                (cluster, rejected, ignored) = handle_pseudogenes(overlapping_pseudo)
            split_clusters, misc = split_cluster(cluster)
            for cluster in split_clusters:
                cluster, rejected, ignored = handle_rfam_only(cluster)
                cluster.validate_members()
                state.reject_all(rejected)
                state.ignore_all(ignored)

                if cluster.rna_type.is_a("rRNA"):
                    clusters = rrna.classify_cluster(cluster)
                    state.add_clusters(cluster)
                else:
                    state.ignore_all(cluster.locations())
        yield state.finalize()


def from_json(handle) -> ty.Iterable[data.FinalizedState]:
    clusters = cluster(handle)
    return build(clusters, IntervalTree())
