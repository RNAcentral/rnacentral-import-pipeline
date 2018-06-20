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

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from .description import description_of

from .rna_type import rna_type_of


def database_names(accessions):
    """
    Generates a comma separated list of all database names in the given list of
    accesions.
    """

    names = {acc.pretty_database for acc in accessions}
    return ','.join(sorted(names))


def partioned_accessions(all_accessions):
    """
    Parition the list of accessions between those athat are active and those
    that are inactive. This returns a tuple of active, inactive Accession
    objects.
    """

    accessions = []
    inactive_accessions = []
    for accession in all_accessions:
        if accession.pop('is_active'):
            accessions.append(Accession.build(accession))
        else:
            inactive_accessions.append(Accession.build(accession))
    return accessions, inactive_accessions


@attr.s()
class Accession(object):
    """
    This represents the minimal amount of information we need to produce a good
    name from the accession level data for a sequence.
    """

    gene = attr.ib(validator=optional(is_a(basestring)))
    optional_id = attr.ib(validator=optional(is_a(basestring)))
    pretty_database = attr.ib(validator=is_a(basestring))
    feature_name = attr.ib(validator=is_a(basestring))
    ncrna_class = attr.ib(validator=optional(is_a(basestring)))
    species = attr.ib(validator=optional(is_a(basestring)))
    common_name = attr.ib(validator=optional(is_a(basestring)))
    description = attr.ib(validator=is_a(basestring))
    locus_tag = attr.ib(validator=optional(is_a(basestring)))

    @classmethod
    def build(cls, data):
        """
        Create a new Accession from the given dict. This assumes the dict has
        keys with the same names as accession fields.
        """
        return cls(**data)  # pylint: disable=star-args

    @property
    def database(self):
        """
        The normalized (lowercase) database name.
        """
        return self.pretty_database.lower()  # pylint: disable=no-member

    @property
    def rna_type(self):
        """
        Get a single INSDC RNA type for this accession.
        """

        if self.feature_name == 'ncRNA':
            return self.ncrna_class
        return self.feature_name


@attr.s()
class Sequence(object):
    """
    The base class that SpeciesSpecific and GenericSequences inherit from. This
    is the representation of sequences for the purpose of our precompute step.
    A Generic sequence represents the information about a single sequence
    across all taxids that sequence has been observed in, while a species one
    is specific to a taxid.
    """

    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=optional(is_a(int)))
    length = attr.ib(validator=is_a(int))
    accessions = attr.ib(validator=is_a(list))
    inactive_accessions = attr.ib(validator=is_a(list))
    is_active = attr.ib(validator=is_a(bool))
    xref_has_coordinates = attr.ib(validator=is_a(bool))
    rna_was_mapped = attr.ib(validator=is_a(bool))
    previous_data = attr.ib(validator=is_a(dict))

    def is_species_specific(self):
        """
        Check if this sequence is specific to a species or is generic.
        """
        return self.taxid is not None


@attr.s()
class SpeciesSequence(Sequence):
    """
    This represents what data is needed in order to precompute the data for a
    species specific sequence.
    """

    @classmethod
    def build(cls, data):
        """
        Given a dictonary of result from the precompute query this will build a
        SpeciesSequence.
        """

        active, inactive = partioned_accessions(data['accessions'])
        return cls(
            upi=data['upi'],
            taxid=data['taxid'],
            length=data['length'],
            accessions=active,
            inactive_accessions=inactive,
            is_active=any(not d for d in data['deleted']),
            xref_has_coordinates=any(data['xref_has_coordinates']),
            rna_was_mapped=data['rna_was_mapped'],
            previous_data=data['previous'][0],
        )

@attr.s()
class GenericSequence(Sequence):

    @classmethod
    def build(cls, sequences):
        """
        Build a Seqeunce object which represents the generic sequence (no
        assigned taxid) given a list of species specific sequences (have an
        assigend taxid).
        """

        inactive = []
        accessions = []
        for seq in sequences:
            accessions.extend(seq.accessions)
            inactive.extend(seq.inactive_accessions)

        has_coordinates = any(s.xref_has_coordinates for s in sequences)
        return cls(
            upi=sequences[0].upi,
            taxid=None,
            length=sequences[0].length,
            accessions=accessions,
            inactive_accessions=inactive,
            is_active=any(seq.is_active for seq in sequences),
            xref_has_coordinates=has_coordinates,
            rna_was_mapped=any(s.rna_was_mapped for s in sequences),
            previous_data={},
        )


@attr.s()
class Update(object):
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=optional(is_a(int)))
    is_active = attr.ib(validator=is_a(bool))
    rna_type = attr.ib(validator=is_a(basestring))
    description = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(basestring))
    has_coordinates = attr.ib(validator=is_a(bool))
    short_description = attr.ib(validator=is_a(basestring), default='')

    @property
    def rna_id(self):
        if self.taxid is not None:
            return '%s_%i' % (self.upi, self.taxid)
        return self.upi

    def as_writeable(self):
        return [
            self.rna_id,
            self.upi,
            self.taxid,
            self.is_active,
            self.description,
            self.rna_type,
            self.has_coordinates,
            self.databases,
        ]


class ActiveUpdate(Update):
    @classmethod
    def build(cls, sequence):
        rna_type = rna_type_of(sequence)
        description = description_of(rna_type, sequence).encode('utf-8')

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
        )


class InactiveUpdate(Update):
    @classmethod
    def build(cls, sequence):
        has_coordinates = sequence.xref_has_coordinates or \
            sequence.rna_was_mapped

        if 'rna_type' in sequence.previous_data:
            rna_type = sequence.previous_data['rna_type']
        else:
            rna_types = {acc.rna_type for acc in sequence.inactive_accessions}
            if len(rna_types) == 1:
                rna_type = rna_types.pop()
            else:
                rna_type = 'ncRNA'

        if 'description' in sequence.previous_data:
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
        )
