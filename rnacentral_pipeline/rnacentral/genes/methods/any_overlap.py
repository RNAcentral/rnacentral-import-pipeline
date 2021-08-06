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

from rnacentral_pipeline.rnacentral.genes.data import State, LocationInfo, Context
from .common_filter import filter


class AnyOverlapMethod:
    def handle_location(self, state: State, context: Context, location: LocationInfo):
        if filter(state, context, location):
            return

        overlapping = state.overlaps(location)
        if not overlapping:
            state.add_singleton_cluster(location)
        elif len(overlapping) == 1:
            cluster = overlapping[0]
            state.add_to_cluster(location, cluster.id)
        else:
            cluster_id = state.merge_clusters(overlapping)
            state.add_to_cluster(location, cluster_id)
