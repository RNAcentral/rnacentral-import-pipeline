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
import itertools as it

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from .description import description_of

from .rna_type import rna_type_of


@attr.s()
class Accession(object):
    """
    This represents the minimal amount of information we need to produce a good
    name from the accession level data for a sequence.
    """

    gene = attr.ib()
    optional_id = attr.ib()
    pretty_database = attr.ib(validator=is_a(basestring))
    feature_name = attr.ib()
    ncrna_class = attr.ib()
    species = attr.ib()
    common_name = attr.ib()
    description = attr.ib()

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
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=optional(is_a(int)))
    accessions = attr.ib(validator=is_a(list))
    xref_status = attr.ib(validator=is_a(list))
    xref_has_coordinates = attr.ib(validator=is_a(bool))
    rna_was_mapped = attr.ib(validator=is_a(bool))

    @classmethod
    def species_to_generic(cls, sequences):
        """
        Build a Seqeunce object which represents the generic sequence (no
        assigned taxid) given a list of species specific sequences (have an
        assigend taxid).
        """

        def get_flat(name):
            """
            Get a flattend list for the given attribute name of sequence
            objects.
            """
            getter = op.itemgetter(name)
            return list(it.chain.from_iterable(getter(s) for s in sequences))

        has_coordinates = any(s['xref_has_coordinates'] for s in sequences)

        return cls(
            upi=sequences[0]['upi'],
            taxid=None,
            accessions=[Accession.build(a) for a in get_flat('accessions')],
            xref_status=get_flat('deleted'),
            xref_has_coordinates=has_coordinates,
            rna_was_mapped=any(s['rna_was_mapped'] for s in sequences),
        )

    @classmethod
    def build(cls, data):
        return cls(
            upi=data['upi'],
            taxid=data['taxid'],
            accessions=[Accession.build(a) for a in data['accessions']],
            xref_status=data['deleted'],
            xref_has_coordinates=any(data['xref_has_coordinates']),
            rna_was_mapped=data['rna_was_mapped'],
        )

    def is_species_specific(self):
        return self.taxid is not None


@attr.s
class UpdatedData(object):
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=optional(is_a(int)))
    is_active = attr.ib(validator=is_a(bool))
    rna_type = attr.ib(validator=is_a(basestring))
    description = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(basestring))
    has_coordinates = attr.ib(validator=is_a(bool))
    short_description = attr.ib(validator=is_a(basestring), default='')

    @classmethod
    def build(cls, sequence):
        print(sequence)
        rna_type = rna_type_of(sequence)
        description = description_of(rna_type, sequence)

        has_coordinates = sequence.xref_has_coordinates or \
            sequence.rna_was_mapped

        databases = {acc.pretty_database for acc in sequence.accessions}
        databases = ','.join(sorted(databases))
        return cls(
            upi=sequence.upi,
            taxid=sequence.taxid,
            is_active=any(d == 'N' for d in sequence.xref_status),
            rna_type=rna_type,
            description=description,
            databases=databases,
            has_coordinates=has_coordinates,
            short_description='',
        )

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
