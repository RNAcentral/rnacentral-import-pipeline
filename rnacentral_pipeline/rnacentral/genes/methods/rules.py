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

import logging
import typing as ty

from rnacentral_pipeline.rnacentral.genes.data import (
    State,
    LocationInfo,
    Context,
    Cluster,
)
from .common_filter import filter, filter_overlaps

LOGGER = logging.getLogger(__name__)


class RuleMethod:
    def has_compatible_rna_types(
        self, location: LocationInfo, cluster: Cluster
    ) -> bool:
        rna_types = cluster.rna_types()
        if not rna_types:
            raise ValueError("This should be impossible")

        for rna_type in rna_types:
            if rna_type.normalized_term != location.rna_type.normalized_term:
                return False
        return True

    def select_mergable(
        self, location: LocationInfo, clusters: ty.List[Cluster]
    ) -> ty.Optional[ty.List[Cluster]]:
        to_merge = []
        target = location.as_interval()
        for cluster in clusters:
            if not cluster.as_interval().overlaps(target):
                continue
            if not self.has_compatible_rna_types(location, cluster):
                continue
            to_merge.append(cluster)

        if not to_merge:
            return None
        return to_merge

    def handle_location(self, state: State, context: Context, location: LocationInfo):
        if filter(state, context, location):
            return

        overlaps = state.overlaps(location)
        if not overlaps:
            LOGGER.debug("Adding singleton cluster of %s", location.id)
            state.add_singleton_cluster(location)
            return

        if filter_overlaps(state, context, location, overlaps):
            return

        to_merge = self.select_mergable(location, overlaps)
        if to_merge is None:
            LOGGER.debug("Adding singleton cluster of %s", location.id)
            state.add_singleton_cluster(location)
            return
        elif len(to_merge) == 1:
            cluster = to_merge[0]
            LOGGER.debug("Adding location %i to cluster %i", location.id, cluster.id)
            state.add_to_cluster(location, cluster.id)
        else:
            LOGGER.debug("Merging into %s", [c.id for c in to_merge])
            cluster_id = state.merge_clusters(to_merge)
            state.add_to_cluster(location, cluster_id)
