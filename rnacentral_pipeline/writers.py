# -*- coding: utf-8 -*-

from __future__ import annotations

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

import csv
import typing as ty
from pathlib import Path

from contextlib import ExitStack
from contextlib import contextmanager

import attr

from rnacentral_pipeline.databases import data


_C = ty.TypeVar("_C")


@attr.s()
class EntryWriter:
    accessions = attr.ib()
    short_sequences = attr.ib(
        metadata={"csv_options": {"delimiter": ",", "lineterminator": "\n"}},
    )
    long_sequences = attr.ib(
        metadata={"csv_options": {"delimiter": ",", "lineterminator": "\n"}}
    )
    references = attr.ib()
    ref_ids = attr.ib()
    regions = attr.ib()
    secondary_structure = attr.ib()
    related_sequences = attr.ib()
    features = attr.ib()
    interactions = attr.ib()
    terms = attr.ib()
    go_annotations = attr.ib()
    go_publication_mappings = attr.ib()

    def write(self, entries: ty.Iterable[data.Entry]):
        for entry in entries:
            if not entry.is_valid():
                continue

            self.accessions.writerows(entry.write_ac_info())
            self.short_sequences.writerows(entry.write_seq_short())
            self.long_sequences.writerows(entry.write_seq_long())
            self.references.writerows(entry.write_refs())
            self.ref_ids.writerows(entry.write_ref_ids())
            self.regions.writerows(entry.write_sequence_regions())
            self.secondary_structure.writerows(entry.write_secondary_structure())
            self.related_sequences.writerows(entry.write_related_sequences())
            self.features.writerows(entry.write_sequence_features())
            self.interactions.writerows(entry.write_interactions())
            self.terms.writerows(entry.write_ontology_terms())

            for annotations in entry.go_annotations:
                self.go_annotations.writerows(annotations.writeable())
                self.go_publication_mappings.writerows(
                    annotations.writeable_publication_mappings()
                )
                self.terms.writerows(annotations.writeable_ontology_terms())


@attr.s()
class OntologyAnnnotationWriter:
    go_annotations = attr.ib()
    go_publication_mappings = attr.ib()
    terms = attr.ib()

    def write(self, annotations: ty.Iterable[data.GoTermAnnotation]):
        for anno in annotations:
            self.go_annotations.writerows(anno.writeable())
            self.go_publication_mappings.writerows(
                anno.writeable_publication_mappings()
            )
            self.terms.writerows(anno.writeable_ontology_terms())


@contextmanager
def build(cls: ty.Type[_C], path: Path) -> ty.Iterator[_C]:
    handles = {}
    with ExitStack() as stack:
        for field in attr.fields(cls):
            out = path / f"{field.name}.csv"
            handle = Path(out).open("w")
            options = field.metadata.get(
                "csv_options",
                {
                    "delimiter": ",",
                    "quotechar": '"',
                    "quoting": csv.QUOTE_ALL,
                    "lineterminator": "\n",
                },
            )
            handles[field.name] = csv.writer(handle, **options)
            stack.enter_context(handle)
        yield cls(**handles)


def entry_writer(path: Path) -> ty.ContextManager[EntryWriter]:
    return build(EntryWriter, path)
