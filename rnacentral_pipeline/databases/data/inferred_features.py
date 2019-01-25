# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import json
import collections as coll

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from more_itertools import windowed


def exonic_features(region, sequence_length):
    """
    Infer the features that represent the exon/intron junctions for a
    particular region in an entry.
    """

    if len(region.exons) <= 1:
        return

    offset = 0
    zero = region.start
    for (first, second) in windowed(region.exons, 2):
        start = (first.stop - zero) + offset
        stop = start + 1

        if stop > sequence_length:
            raise ValueError("%i is beyond sequence end: %i" % (stop, sequence_length))

        yield InferredSequenceFeature(
            start,
            stop,
            "exon_junction",
        )
        offset += first.length()


@attr.s(frozen=True, hash=True, cmp=True)
class InferredSequenceFeature(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))
    relationship = attr.ib(validator=is_a(basestring))
    metadata = attr.ib(
        validator=optional(is_a(dict)),
        factory=dict,
        hash=False,
        cmp=False,
    )

    def writeable(self, *prefix):
        """
        Generate an array to write out to represent this
        InferredSequenceFeature object. This is mean to be written to a file
        and the loaded into the database.
        """

        data = list(prefix)
        data.extend([
            self.start,
            self.stop,
            self.relationship,
            json.dumps(self.metadata)
        ])
        return data


class EntryFeatureInference(object):
    """
    A class to infer sequence features for entries.
    """

    def infer_related_features(self, entry):
        """
        Determine the features based off the related sequences.
        """

        for related in entry.related_sequences:
            for endpoints in related.coordinates:
                metadata = {'related': related.sequence_id}
                yield InferredSequenceFeature(
                    endpoints.start,
                    endpoints.stop,
                    related.relationship,
                    metadata,
                )

    def infer_exon_junctions(self, entry):
        """
        Infer all Exon/Intron junction features. This will examine the
        SequenceRegions to determine the correct features. Each feature will
        represent a single exon/intron junction which occurs in one or more
        region. If there are several regions with the same junctions only one
        set of features will be created.
        """

        junctions = coll.defaultdict(set)
        for region in entry.regions:
            for feature in exonic_features(region, entry.sequence_length):
                junctions[feature].add(region.name)

        for junction, _ in junctions.items():
            if junction.start == 0 and junction.stop == len(entry.sequence):
                continue

            yield attr.evolve(
                junction,
                metadata={},  # No region tracking yet
            )

    def features(self, entry):
        """
        Infer all features for the given entry. This will return a iterable of
        InferredSequenceFeature objects which represents the features for the
        entry.
        """

        methods = [
            self.infer_related_features,
            self.infer_exon_junctions,
        ]

        for method in methods:
            for feature in method(entry):
                yield feature

    def writeables(self, entry):
        """
        Generate an iterable of all writeable arrays of inferred features for
        the given entry.
        """
        for feature in self.features(entry):
            yield feature.writeable(entry.accession, entry.ncbi_tax_id)


class HitFeatureInference(object):
    """
    Compute sequence features based off the exon/intron structure of a sequence
    match.
    """

    def features(self, hit):
        """
        Infer all features for the given entry. This will return a iterable of
        InferredSequenceFeature objects which represents the features for the
        entry.
        """

        for feature in exonic_features(hit, hit.sequence_length):
            yield feature

    def writeables(self, hit):
        """
        Generate an iterable of all writeable arrays of inferred features for
        the given entry.
        """

        for feature in self.features(hit):
            yield feature.writeable(hit.upi)
