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

import operator as op
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline.databases.data.utils import INSDC_SO_MAPPING
from rnacentral_pipeline.rnacentral.precompute import qa
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.description import (
    description_of, short_description_for)
from rnacentral_pipeline.rnacentral.precompute.rna_type import rna_type_of
from rnacentral_pipeline.rnacentral.repeats import tree


@attr.s(frozen=True)
class SequenceUpdate:
    """
    This represents the data that is an update to our precomputed data.
    """

    sequence = attr.ib(validator=is_a(Sequence))
    insdc_rna_type = attr.ib(validator=is_a(str))
    so_rna_type = attr.ib(validator=is_a(str))
    description = attr.ib(validator=is_a(str))
    short_description = attr.ib(validator=is_a(str))
    qa_status = attr.ib(validator=is_a(qa.QaStatus))

    @classmethod
    def active(cls, context: Context, sequence: Sequence) -> "SequenceUpdate":
        insdc_rna_type = rna_type_of(sequence)
        so_rna_type = INSDC_SO_MAPPING.get(insdc_rna_type, "SO:0000655")
        description = description_of(insdc_rna_type, sequence)
        short_description = short_description_for(description, sequence)

        return cls(
            sequence=sequence,
            insdc_rna_type=insdc_rna_type,
            so_rna_type=so_rna_type,
            description=description,
            short_description=short_description,
            qa_status=qa.status(context, sequence, insdc_rna_type),
        )

    @classmethod
    def inactive(cls, sequence: Sequence) -> "SequenceUpdate":
        """
        This will build InactiveUpdate for the given sequence. This will try to
        copy over any previous data and will do very little work to create a
        correct RNA type or useful description otherwise.
        """

        insdc_rna_type = sequence.previous_update.get("rna_type", "")
        if not insdc_rna_type:
            insdc_rna_type = "ncRNA"
            insdc_rna_types: ty.Set[str] = {
                acc.rna_type for acc in sequence.inactive_accessions
            }
            if len(insdc_rna_types) == 1:
                insdc_rna_type = insdc_rna_types.pop()

        description = sequence.previous_update.get("description", "")
        if not description:
            description = "{name} {insdc_rna_type}"

        so_rna_type = INSDC_SO_MAPPING.get(insdc_rna_type, "SO:0000655")
        return cls(
            sequence=sequence,
            insdc_rna_type=insdc_rna_type,
            so_rna_type=so_rna_type,
            description=description,
            short_description=description,
            qa_status=None,
        )

    @classmethod
    def from_sequence(cls, context: Context, sequence: Sequence) -> "SequenceUpdate":
        if sequence.is_active:
            return cls.active(context, sequence)
        return cls.inactive(sequence)

    @property
    def is_active(self):
        return self.sequence.is_active

    def database_names(self) -> str:
        """
        Generates a comma separated list of all database names in the given list of
        accesions.
        """

        databases = {acc.pretty_database for acc in self.sequence.accessions}
        return ",".join(sorted(databases, key=op.methodcaller("lower")))

    def as_writeables(self) -> ty.Iterable[ty.List[str]]:
        """
        Yield the arrays that will be written out for updates.
        """

        yield [
            self.sequence.rna_id,
            self.sequence.upi,
            str(self.sequence.taxid),
            str(int(self.sequence.is_active)),
            self.description,
            self.insdc_rna_type,
            str(int(self.sequence.has_coordinates)),
            self.database_names(),
            self.short_description,
            str(self.sequence.last_release),
            self.so_rna_type,
        ]

    def writeable_statuses(self) -> ty.Iterable[ty.List[str]]:
        """
        Yield arrays for all updates for qa status.
        """

        yield self.qa_status.writeable(self.sequence.upi, self.sequence.taxid)


@attr.s(frozen=True)
class GenericUpdate:
    upi = attr.ib(validator=is_a(str))
    updates: ty.List[SequenceUpdate] = attr.ib(validator=is_a(list))
    is_active = attr.ib(validator=is_a(bool))

    @classmethod
    def active(cls, updates: ty.List[SequenceUpdate]) -> "GenericUpdate":
        return cls(upi=updates[0].sequence.upi, updates=updates, is_active=True,)

    @classmethod
    def inactive(cls, updates: ty.List[SequenceUpdate]) -> "GenericUpdate":
        return cls(upi=updates[0].sequence.upi, updates=updates, is_active=False,)

    @classmethod
    def from_updates(
        cls, context: Context, updates: ty.List[SequenceUpdate]
    ) -> "GenericUpdate":
        active_updates = [u for u in updates if u.is_active]
        if active_updates:
            return cls.active(active_updates)
        return cls.inactive(updates)

    @property
    def last_release(self):
        return max(u.sequence.last_release for u in self.updates)

    @property
    def species_count(self):
        all_species = set()
        for update in self.updates:
            all_species.update(update.sequences.species())
        return len(all_species)

    @property
    def description(self) -> str:
        return f"{self.insdc_rna_type} from {self.species_count} species"

    @property
    def database_names(self) -> str:
        databases: ty.Set[str] = set()
        for update in self.updates:
            for accession in update.sequence.accessions:
                databases.add(accession.pretty_database)
        return ",".join(sorted(databases, key=op.methodcaller("lower")))

    @property
    def insdc_rna_type(self) -> str:
        rna_types = {u.insdc_rna_type for u in self.updates}
        if len(rna_types) == 1:
            return rna_types.pop()
        return "ncRNA"

    @property
    def so_rna_type(self) -> str:
        rna_types = {u.so_rna_type for u in self.updates}
        if len(rna_types) == 1:
            return rna_types.pop()
        return "SO:0000655"

    def as_writeables(self) -> ty.Iterable[ty.List[str]]:
        yield [
            self.upi,
            self.upi,
            "",
            str(int(self.is_active)),
            self.description,
            self.insdc_rna_type,
            "",
            self.database_names,
            self.description,
            str(self.last_release),
            self.so_rna_type,
        ]

    def writeable_status(self) -> ty.Iterable[ty.List[str]]:
        pass
