# -*- coding: utf-8 -*-

from __future__ import annotations

"""
Copyright [2009-${2022}] EMBL-European Bioinformatics Institute
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

import enum
import typing as ty
from pathlib import Path
import operator as op
import itertools as it

import attr
from pypika import Query, Table
from pypika import functions as fn
import pandas as pd
# import sklearn as skl
import networkx as nx

from rnacentral_pipeline.rnacentral.genes import features


@attr.s
class Exon:
    start: int
    stop: int


@attr.s
class Region:
    name: str
    group: str
    exons: ty.List[Exon]
    length: int
    so_rna_type: str

    @property
    def start(self):
        return min(e.start for e in self.exons)

    @property
    def stop(self):
        return min(e.stop for e in self.exons)

    def exon_endpoints(self):
        return sorted([(e.start, e.stop) for e in self.exons])


# Encode as both SO types
# Encode similarity between SO types
# Encode simiarity between exon/start stop with various measures
# Try decision tree vs random forest vs Support vector machine

@attr.frozen
class DataPoint:
    region1: Region
    region2: Region
    should_join: bool

    @classmethod
    def build_joined(cls, left: Region, right: Region) -> DataPoint:
        return cls(
            region1=left,
            region2=right,
            should_join=True
        )

    @classmethod
    def build_distinct(cls, left: Region, right: Region) -> DataPoint:
        return cls(
            region1=left,
            region2=right,
            should_join=False
        )

class Metadata:
    so_tree: nx.Graph
    allowed_so_terms: ty.Set[str]


@enum.unique
class Feature(enum.Enum):
    REGION_JACCARD = "region_jaccard_index"
    EXON_JUNC_JACCARD = "exon_junction_jaccard_index"
    LENGTH = "length"
    LEFT_SO_RNA_TYPE = "left_so_rna_type"
    RIGHT_SO_RNA_TYPE = "right_so_rna_type"
    SO_TERM_TREE_DISTANCE = "so_term_tree_distance"

    def compute(self, metadata: Metadata, left: Region, right: Region):
        if self is Feature.REGION_JACCARD:
            return features.region_jaccard(left, right)
        if self is Feature.EXON_JUNC_JACCARD:
            return features.exon_junction_jaccard_index(left, right)
        if self is Feature.LEFT_SO_RNA_TYPE:
            return left.so_rna_type
        if self is Feature.RIGHT_SO_RNA_TYPE:
            return right.so_rna_type
        if self is Feature.SO_TERM_TREE_DISTANCE:
            assert left.so_rna_type in metadata.allowed_so_terms
            assert right.so_rna_type in metadata.allowed_so_terms
            return features.tree_distance(metadata.so_tree, left.so_rna_type, right.so_rna_type)
        raise ValueError("Not yet implemented")


@attr.s
class Parameters:
    test_train_split: float
    distance_to_negative_example: int
    features: ty.List[ty.Tuple[str, Feature]]

    @classmethod
    def default(cls) -> Parameters:
        return cls(
            test_train_split=0.3,
            distance_to_negative_example=10000,
            features=[
                ('exon_distance', Feature.EXON_JUNC_JACCARD),
                ('so_measure', Feature.SO_TERM_TREE_DISTANCE),
            ],
        )

REGIONS = Table('rnc_sequence_regions')
EXONS = Table('rnc_sequence_region_exons')
XREF = Table('xref')
ACC = Table('rnc_accessions')


def query_by_gene(genes) -> ty.Iterable[Query]:
    yield Query.from_(REGIONS).\
            select(REGIONS.urs_taxid, ACC.gene, REGIONS.length, EXONS.start, EXONS.stop).\
            join(EXONS, EXONS.region_id == REGIONS.id).\
            join(XREF, fn.Concat(XREF.urs, "_", XREF.taxid) == REGIONS.urs_taxid).\
            join(ACC, ACC.accession == XREF.ac).\
            where(
                ACC.gene.in_(genes)
            )


def queries() -> ty.Iterable[Query]:
    yield from query_by_gene()


def overlaps(left: Region, right: Region) -> bool:
    return left.start <= right.start and left.stop >= right.stop


def close_enough(left: Region, right: Region, distance=1000) -> bool:
    return abs(left.start - right.start) < distance or abs(left.stop - right.stop) < distance


def trainable_data(data: ty.Iterable[Region], parameters: Parameters=Parameters.default()) -> ty.Iterable[DataPoint]:
    grouped = it.groupby(data, op.attrgetter('group'))
    grouped = {g: list(v) for g, v in grouped}
    groups = sorted(grouped.keys())
    for group1, group2 in it.product(groups, groups):
        for (left, right) in it.product(grouped[group1], grouped[group2]):
            if left == right:
                continue
            if group1 == group2:
                yield DataPoint.build_joined(left, right)
            elif overlaps(left, right) or close_enough(left, right, distance=parameters.distance_to_negative_example):
                yield DataPoint.build_distinct(left, right)


def training_data(metadata: Metadata, points: ty.Iterable[DataPoint], params: Parameters=Parameters.default()) -> pd.DataFrame:
    data = []
    for point in points:
        entry = {}
        for (name, feature) in params.features:
            entry[name] = feature.compute(metadata, point.region1, point.region2)
        data.append(entry)
    return pd.DataFrame(data)


# def train_decision_tree(training: pd.DataFrame) -> skl.tree.DecisionTreeClassifier:
#     clf = skl.tree.DecisionTreeClassifier()
#     clf.fit(training.drop("response", axis=1), training["response"])
#     return clf
