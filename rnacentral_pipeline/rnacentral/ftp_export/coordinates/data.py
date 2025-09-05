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

import itertools as it
import json
import operator as op
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline.databases.data import Database, regions


def lookup_databases(raw):
    databases = []
    for db in raw:
        if db is None:
            continue
        found = Database.lookup(db)
        if not found:
            raise ValueError(f"Failed to lookup database {db}")
        databases.append(found.pretty())
    return databases


def clean_databases(raw):
    return [d.replace(" ", "_") for d in raw]


@attr.s(hash=True, slots=True, frozen=True)
class Region(object):
    region_id = attr.ib(validator=is_a(str), converter=str)
    rna_id = attr.ib(validator=is_a(str), converter=str)
    region = attr.ib(validator=is_a(regions.SequenceRegion))
    was_mapped = attr.ib(validator=is_a(bool))
    is_gene = attr.ib(validator=is_a(bool))
    gene_start = attr.ib(validator=optional(is_a(int)), default=None, cmp=False)
    gene_end = attr.ib(validator=optional(is_a(int)), default=None, cmp=False)
    identity = attr.ib(
        validator=optional(is_a(float)),
        default=None,
        cmp=False,
    )
    metadata = attr.ib(validator=is_a(dict), default=dict, cmp=False)

    @classmethod
    def build(cls, index, raw):
        is_gene = False
        identity = None
        exons = []
        metadata = {}
        rna_id = None
        gene_start = None
        gene_end = None

        if raw["rna_id"] is None:
            is_gene = True

        ## Treat genes and transcripts separately
        if is_gene:
            region_id = raw["gene_name"]

            metadata["description"] = raw["description"]
            metadata["rna_type"] = raw["rna_type"]
            metadata["providing_databases"] = []
            metadata["databases"] = raw["databases"]

            gene_start = raw["gene_start"]
            gene_end = raw["gene_end"]

            rna_id = raw["gene_name"]
        else:
            if raw["identity"] is not None:
                identity = float(raw["identity"])

            for exon in raw["exons"]:
                exons.append(
                    regions.Exon(
                        start=exon["exon_start"],
                        stop=exon["exon_stop"],
                    )
                )

            region_id = "{rna_id}.{index}".format(rna_id=raw["rna_id"], index=index)

            metadata["description"] = raw["description"]
            metadata["rna_type"] = raw["rna_type"]
            metadata["providing_databases"] = lookup_databases(
                raw["providing_databases"]
            )
            metadata["databases"] = clean_databases(raw["databases"])
            rna_id = raw["rna_id"]

            if raw["gene_name"] is not None:
                metadata["parent_gene"] = raw["gene_name"]

        ## End separate treatment, should be able to be unified now

        if not metadata["providing_databases"]:
            if not raw["was_mapped"] and not is_gene:
                raise ValueError(
                    f"No providing database for an unmapped region!\n{raw}"
                )
            del metadata["providing_databases"]

        return cls(
            region_id=region_id,
            rna_id=rna_id,
            region=regions.SequenceRegion(
                assembly_id=raw["assembly_id"],
                chromosome=raw["chromosome"],
                strand=raw["strand"],
                exons=exons,
                coordinate_system=regions.CoordinateSystem.one_based(),
            ),
            identity=identity,
            was_mapped=raw["was_mapped"],
            is_gene=is_gene,
            metadata=metadata,
            gene_start=gene_start,
            gene_end=gene_end,
        )

    @property
    def start(self):
        if self.is_gene:
            return self.gene_start
        return self.region.start

    @property
    def stop(self):
        if self.is_gene:
            return self.gene_end
        return self.region.stop

    @property
    def chromosome(self):
        return self.region.chromosome

    @property
    def exons(self):
        return self.region.exons

    @property
    def source(self):
        if self.was_mapped:
            return "alignment"
        if self.is_gene:
            return "gene-prediction"
        return "expert-database"

    def string_strand(self):
        return self.region.strand.display_string()

    def as_one_based(self):
        return attr.evolve(self, region=self.region.as_one_based())

    def as_zero_based(self):
        return attr.evolve(self, region=self.region.as_zero_based())


def parse(regions) -> ty.Iterable[Region]:
    for index, region in enumerate(regions):
        yield Region.build(index, region)


def from_file(handle) -> ty.Iterable[Region]:
    data = map(json.loads, handle)
    data = filter(lambda d: d["rna_type"] != "NULL", data)
    yield from parse(data)
