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

import operator as op

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from .. import qa
from ..description import description_of
from ..description import short_description_for
from ..rna_type import rna_type_of


def database_names(accessions):
    """
    Generates a comma separated list of all database names in the given list of
    accesions.
    """

    databases = {acc.pretty_database for acc in accessions}
    return ','.join(sorted(databases, key=op.methodcaller('lower')))


@attr.s()
class QaStatus(object):
    """
    This represents an update to the QA status table.
    """

    incomplete_sequence = attr.ib(validator=is_a(bool))
    possible_contamination = attr.ib(validator=is_a(bool))
    missing_rfam_match = attr.ib(validator=is_a(bool))
    is_repetitive = attr.ib(validator=is_a(bool))
    no_data = attr.ib(validator=optional(is_a(bool)), default=False)

    @classmethod
    def build(cls, rna_type, data):
        """
        This will build a new QaStatus object update given the particular RNA
        type and a Sequence object.
        """
        return cls(
            incomplete_sequence=qa.incomplete_sequence(rna_type, data),
            possible_contamination=qa.possible_contamination(rna_type, data),
            missing_rfam_match=qa.missing_rfam_match(rna_type, data),
            is_repetitive=False,
        )

    @classmethod
    def empty(cls):
        """
        Build a QaStatus that will not produce any writeable data. Objects
        produced by this should not be used for any updates.
        """
        return cls(False, False, False, False, no_data=True)

    @property
    def has_issue(self):
        """
        Check if this QA update indicates if there is any issue.
        """

        fields = attr.fields(self.__class__)
        return any(getattr(self, field.name) for field in fields)

    def writeable(self, upi, taxid):
        """
        Create a writeable array for writing CSV files.
        """

        if self.no_data:
            return None

        return [
            '%s_%i' % (upi, taxid),
            upi,
            taxid,
            int(self.has_issue),
            int(self.incomplete_sequence),
            int(self.possible_contamination),
            int(self.missing_rfam_match),
        ]


@attr.s()
class Update(object):
    """
    This represents the data that is an update to our precomputed data.
    """

    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=optional(is_a(int)))
    is_active = attr.ib(validator=is_a(bool))
    rna_type = attr.ib(validator=is_a(basestring))
    description = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(basestring))
    has_coordinates = attr.ib(validator=is_a(bool))
    qa_status = attr.ib(validator=is_a(QaStatus))
    short_description = attr.ib(validator=is_a(basestring))

    @property
    def rna_id(self):
        """
        Build the RNA id which is {upi}_{taxid}.
        """

        if self.taxid is not None:
            return '%s_%i' % (self.upi, self.taxid)
        return self.upi

    def as_writeables(self):
        """
        Yield the arrays that will be written out for updates.
        """
        yield [
            self.rna_id,
            self.upi,
            self.taxid,
            int(self.is_active),
            self.description,
            self.rna_type,
            int(self.has_coordinates),
            self.databases,
        ]

    def writeable_statuses(self):
        """
        Yield arrays for all updates for qa status.
        """

        # pylint: disable=no-member
        if self.taxid:
            writeable = self.qa_status.writeable(self.upi, self.taxid)
            if writeable:
                yield writeable


class ActiveUpdate(Update):
    """
    This is an update for active sequences.
    """

    @classmethod
    def build(cls, sequence):
        """
        Build an ActiveUpdate object for the given sequence. This will compute
        the rna_type, description and such data.
        """

        rna_type = rna_type_of(sequence)
        description = description_of(rna_type, sequence).encode('utf-8')
        short_description = short_description_for(description, sequence)

        has_coordinates = sequence.xref_has_coordinates or \
            sequence.rna_was_mapped

        return cls(
            upi=sequence.upi,
            taxid=sequence.taxid,
            is_active=True,
            rna_type=rna_type,
            description=description,
            databases=database_names(sequence.accessions),
            has_coordinates=has_coordinates,
            qa_status=QaStatus.build(rna_type, sequence),
            short_description=short_description,
        )


class InactiveUpdate(Update):
    """
    This is an update for inactive sequences.
    """

    @classmethod
    def build(cls, sequence):
        """
        This will build InactiveUpdate for the given sequence. This will try to
        copy over any previous data and will do very little work to create a
        correct RNA type or useful description otherwise.
        """

        has_coordinates = sequence.xref_has_coordinates or \
            sequence.rna_was_mapped

        if sequence.previous_data and 'rna_type' in sequence.previous_data:
            rna_type = sequence.previous_data['rna_type']
        else:
            rna_types = {acc.rna_type for acc in sequence.inactive_accessions}
            if len(rna_types) == 1:
                rna_type = rna_types.pop()
            else:
                rna_type = 'ncRNA'

        if sequence.previous_data and 'description' in sequence.previous_data:
            description = sequence.previous_data['description']
        else:
            description = '{rna_type} from {count} species'.format(
                rna_type=rna_type,
                count=len({ac.species for ac in sequence.inactive_accessions})
            )

        return cls(
            upi=sequence.upi,
            taxid=sequence.taxid,
            is_active=False,
            rna_type=rna_type,
            description=description,
            databases=database_names(sequence.inactive_accessions),
            has_coordinates=has_coordinates,
            qa_status=QaStatus.empty(),
            short_description='',
        )
