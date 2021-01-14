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

import enum
import logging
import typing as ty

import attr
import networkx as nx
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data import Database, RnaType
from rnacentral_pipeline.databases.sequence_ontology import tree
from rnacentral_pipeline.rnacentral.precompute.data import context
from rnacentral_pipeline.rnacentral.precompute.data import sequence as seq
from rnacentral_pipeline.rnacentral.precompute.qa import contamination as cont

LOGGER = logging.getLogger(__name__)

ACCEPTED_DATABASES = {
    Database.five_srrnadb,
    Database.flybase,
    Database.gtrnadb,
    Database.lncbase,
    Database.lncipedia,
    Database.mirbase,
    Database.mirgenedb,
    Database.pirbase,
    Database.pombase,
    # Database.sgd, This is excluded because they tend ot use things like rRNA_gene, which we don't like as much as the transcript terms
    Database.snodb,
    Database.snorna_database,
    Database.tarbase,
    Database.zwd,
    Database.pdbe,
    Database.gencode,
}

MODEL_GROUPS = {
    "SO:0000209": {
        "RF02543",
        "RF00002",
        "RF01960",
        "RF02546",
        "RF02543",
        "RF02541",
        "RF02540",
        "RF01959",
        "RF00177",
        "RF01960",
        "RF02542",
        "RF02542",
        "RF00957",
    },
}


@enum.unique
class SourceName(enum.Enum):
    generic_database = enum.auto()
    mod = enum.auto()
    r2dt_hit = enum.auto()
    rfam_hit = enum.auto()

    def is_hit(self):
        if self is SourceName.generic_database:
            return False
        if self is SourceName.mod:
            return False
        if self is SourceName.r2dt_hit:
            return True
        if self is SourceName.rfam_hit:
            return True
        raise ValueError(f"Unhandled SourceName {self}")

    def is_database(self):
        if self is SourceName.generic_database:
            return True
        if self is SourceName.mod:
            return True
        if self is SourceName.r2dt_hit:
            return False
        if self is SourceName.rfam_hit:
            return False
        raise ValueError(f"Unhandled SourceName {self}")


@attr.s(hash=True, frozen=True)
class Source:
    name = attr.ib(validator=is_a(SourceName))
    data = attr.ib()

    @classmethod
    def from_mod(cls, accession):
        return cls(name=SourceName.mod, data=accession)

    @classmethod
    def from_generic(cls, accession):
        return cls(name=SourceName.generic_database, data=accession)

    @classmethod
    def from_r2dt(cls, hit):
        return cls(name=SourceName.r2dt_hit, data=hit)

    @classmethod
    def from_rfam(cls, hit):
        return cls(name=SourceName.rfam_hit, data=hit)

    def is_hit(self):
        return self.name.is_hit()

    def is_database(self):
        return self.name.is_database()


@attr.s()
class RnaTypeAnnotation:
    path: ty.List[RnaType] = attr.ib(validator=is_a(list))
    source: ty.Set[Source] = attr.ib(validator=is_a(set))

    @classmethod
    def build(
        cls, context: context.Context, rna_type: RnaType, source: Source
    ) -> "RnaTypeAnnotation":
        path = tree.rna_type_tree(context.so_tree, rna_type.so_id)
        path = [RnaType.from_so_term(context.so_tree, so_id) for (so_id, _) in path]
        if path[-1] != rna_type:
            path.append(rna_type)
        return cls(path=path, source={source})

    @classmethod
    def from_accession(cls, context, accession) -> "RnaTypeAnnotation":
        if accession.database in ACCEPTED_DATABASES:
            return cls.build(context, accession.rna_type, Source.from_mod(accession))
        return cls.build(context, accession.rna_type, Source.from_generic(accession))

    @classmethod
    def from_r2dt(cls, context, hit) -> "RnaTypeAnnotation":
        return cls.build(context, hit.model_rna_type, Source.from_r2dt(hit))

    @classmethod
    def from_rfam(cls, context: context.Context, hit) -> "RnaTypeAnnotation":
        return cls.build(context, hit.model_rna_type, Source.from_rfam(hit))

    def equivalent(self, other: "RnaTypeAnnotation") -> bool:
        return self.path == other.path

    def is_parent_of(self, other: "RnaTypeAnnotation") -> bool:
        if len(self.path) >= len(other.path):
            return False
        for (left, right) in zip(self.path, other.path):
            if left != right:
                return False
        return True

    def is_child_of(self, other: "RnaTypeAnnotation") -> bool:
        return other.is_parent_of(self)

    def is_mergeable(self, other: "RnaTypeAnnotation"):
        return (
            self.path == other.path
            or self.is_parent_of(other)
            or self.is_child_of(other)
        )

    @property
    def rna_type(self):
        return self.path[-1]

    def merge(self, other: "RnaTypeAnnotation"):
        source = set()
        source.update(self.source)
        source.update(other.source)
        if self.path == other.path:
            self.source = source
        elif self.is_parent_of(other):
            self.source = source
            self.path = list(other.path)
        elif self.is_child_of(other):
            self.source = source
        else:
            raise ValueError(f"Cannot merge {self.path} to {other.path}")

    def source_matches(self, fn):
        return any(fn(s) for s in self.source)

    def has_source(self, source: SourceName):
        return self.source_matches(lambda s: s.name == source)

    def is_specific(self):
        return len(self.path) > 1


def merge_annotations(given: ty.List[RnaTypeAnnotation]) -> ty.List[RnaTypeAnnotation]:
    grouped = [given[0]]
    for annotation in given[1:]:
        for group in grouped:
            if group.is_mergeable(annotation):
                group.merge(annotation)
                break
        else:
            grouped.append(annotation)
    return grouped


def covers_sequence(hit) -> bool:
    return hit.sequence_info.completeness > 0.80 and hit.model_info.completeness > 0.80


def is_mito_hit(context, sequence, hit) -> bool:
    return (
        sequence.is_mitochondrial()
        and context.term_is_a("rRNA", hit.model_rna_type)
        and hit.model in cont.ALLOWED_FAMILIES
    )


def mito_rna_type_of(context, _hit) -> RnaType:
    return RnaType.from_so_term(context.so_tree, "SO:0002128")


def is_lncrna_data(
    ctx: context.Context, annotations: ty.List[RnaTypeAnnotation]
) -> bool:
    def allowed(term):
        return ctx.term_is_a("lnc_RNA", term) or ctx.term_is_a("antisense_RNA", term)

    selected = filter(lambda a: a.source_matches(lambda s: s.is_database()), annotations)
    return all(allowed(a.rna_type) for a in selected)


def has_similar_rfam_and_r2dt(context, annotations):
    r2dt = [a for a in annotations if a.has_source(SourceName.r2dt_hit)]
    rfam = [a for a in annotations if a.has_source(SourceName.rfam_hit)]
    if len(r2dt) != 1 or len(rfam) != 1:
        return False
    return context.term_is_a('rRNA', r2dt[0].rna_type) and \
        context.term_is_a('rRNA', rfam[0].rna_type)


def rfam_db_annotations(annotations: ty.List[RnaTypeAnnotation]) -> ty.Optional[RnaTypeAnnotation]:
    def from_rfam(annotation: RnaTypeAnnotation) -> bool:
        def fn(source) -> bool:
            return source.name == SourceName.generic_database and \
                    source.data.database == Database.rfam

        return annotation.source_matches(fn)

    selected = [a for a in annotations if from_rfam(a)]
    if len(selected) == 1:
        return selected[0]
    return None


def all_annotations(context: context.Context, sequence: seq.Sequence) -> ty.List[RnaTypeAnnotation]:
    annotations = []
    for accession in sequence.accessions:
        annotations.append(RnaTypeAnnotation.from_accession(context, accession))
    for r2dt in sequence.r2dt_hits:
        if r2dt.paired_ratio() is None or r2dt.paired_ratio() > 0.80:
            annotations.append(RnaTypeAnnotation.from_r2dt(context, r2dt))

    for rfam in sequence.rfam_hits:
        if not covers_sequence(rfam):
            continue
        if is_mito_hit(context, sequence, rfam):
            rfam = attr.evolve(rfam, model_rna_type=mito_rna_type_of(context, rfam))
        ann = RnaTypeAnnotation.from_rfam(context, rfam)
        annotations.append(ann)
    return annotations


def rna_type_of(
    context: context.Context, sequence: seq.Sequence
) -> ty.Optional[RnaType]:

    annotations = all_annotations(context, sequence)
    merged = merge_annotations(annotations)

    if len(merged) == 1:
        return merged[0].rna_type

    mod_annotations = [m for m in merged if m.has_source(SourceName.mod)]
    if len(mod_annotations) == 1 and mod_annotations[0].is_specific():
        return mod_annotations[0].rna_type

    hits = [m for m in merged if m.source_matches(lambda s: s.is_hit())]
    if len(hits) == 1:
        return hits[0].rna_type

    if mod_annotations:
        return mod_annotations[0].rna_type

    if is_lncrna_data(context, merged):
        return RnaType.from_so_id(context.so_tree, "SO:0001877")

    if has_similar_rfam_and_r2dt(context, merged):
        r2dt = [a for a in merged if a.has_source(SourceName.r2dt_hit)]
        if len(r2dt) == 1:
            return r2dt[0].rna_type

    rfam_annotation = rfam_db_annotations(merged)
    if rfam_annotation:
        return rfam_annotation.rna_type

    return RnaType.ncRNA()
