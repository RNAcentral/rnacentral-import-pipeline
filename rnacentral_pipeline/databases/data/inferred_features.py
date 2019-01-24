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

    def writeable(self, entry):
        """
        Generate an array to write out to represent this
        InferredSequenceFeature object. This is mean to be written to a file
        and the loaded into the database.
        """

        return [
            entry.accession,
            entry.ncbi_tax_id,
            self.start,
            self.stop,
            self.relationship,
            json.dumps(self.metadata)
        ]


class FeatureInference(object):
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

    def region_features(self, entry, region):
        """
        Infer the features that represent the exon/intron junctions for a
        particular region in an entry.
        """

        offset = 0
        zero = region.start
        for exon in region.exons:
            start = (exon.start - zero) + offset
            yield InferredSequenceFeature(
                start,
                start + exon.length,
                "exon_junction",
            )
            offset += exon.length

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
            for feature in self.region_features(entry, region):
                junctions[feature].add(region.name)

        for junction, region_ids in junctions.items():
            yield attr.evolve(
                junction,
                metadata={'regions': sorted(region_ids)},
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
            yield feature.writeable(entry)
