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

from rnacentral_pipeline.databases.data.databases import Database
from rnacentral_pipeline.rnacentral.genes.data import State, LocationInfo, Context, Cluster

LOGGER = logging.getLogger(__name__)

def overlaps_pseudogene(context: Context, location: LocationInfo) -> bool:
    return bool(context.pseudogenes.overlap(location.as_interval()))


def overlaps_repetitive_region(context: Context, location: LocationInfo) -> bool:
    return bool(context.repetitive.overlap(location.as_interval()))


# def possible_rfam_overlap(rfam: LocationInfo, other: LocationInfo) -> bool:


def is_rfam_shift(state: State, context: Context, location: LocationInfo, overlaps: ty.List[Cluster]) -> bool:
    return False
#     if not location.is_rfam_only():
#         return False
#     assert not location.has_introns()
#     for cluster in overlaps:
#         locations = state.members_of(cluster.id)
#         rfam_range = rfam.exons[0]
#         for other in locations:
#             if Database.rfam in other.databases:
#                 continue
#             if other.has_introns():
#                 continue
#             other_range = other.exons[0]
#             if (rfam_range.start other_range.start &&
#     return False


def always_bad_location(location: LocationInfo) -> bool:
    if location.qa.has_issue:
        if location.extent.chromosome == "MT":
            state = location.qa.as_tuple()
            # If the only issue is contamination, it is actually good as we
            # misclassified mito sequences sometimes.
            LOGGER.debug("State: %s", state)
            return state != (True, False, True, False)
        return True
    if location.counts.mapped_count > 30:
        return True
    return False


def always_ignorable_location(location: LocationInfo) -> bool:
    return False


def filter(state: State, context: Context, location: LocationInfo) -> bool:
    if always_ignorable_location(location):
        LOGGER.debug("Always Ignoring: %s", location.id)
        state.ignore_location(location)
        return True

    if always_bad_location(location):
        LOGGER.debug("Always Rejected: %s", location.id)
        state.reject_location(location)
        return True

    if overlaps_pseudogene(context, location):
        LOGGER.debug("Rejecting Pseudogene: %s", location.id)
        state.reject_location(location)
        return True

    if overlaps_repetitive_region(context, location):
        LOGGER.debug("Rejected repetitive region: %s", location.id)
        state.reject_location(location)
        return True
    return False


def filter_overlaps(state: State, context: Context, location: LocationInfo, overlaps: ty.List[Cluster]):
    if is_rfam_shift(state, context, location, overlaps):
        LOGGER.debug("Rejecting %s as it is an Rfam shift", location)
        state.reject_location(location)
        return True
    return False
