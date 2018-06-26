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
    organelle = attr.ib(validator=optional(is_a(basestring)))
    lineage = attr.ib(validator=optional(is_a(basestring)))

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

    @property
    def domain(self):
        """
        Get the domain, if any, that is assigned to this accession. This will
        not include uncultured, environmental or synthetic domains.
        """

        if 'uncultured' in self.lineage or \
                'environmental' in self.lineage or \
                'synthetic' in self.lineage:
            return None

        return self.lineage.split(';')[0]  # pylint: disable=no-member

    def is_mitochondrial(self):
        """
        Check if this accession is mitochrondrial.
        """
        return 'mitochondri' in self.description or \
            (self.organelle and 'mitochondri' in self.organelle)



@attr.s()
class RfamHit(object):
    """
    This represents the information needed to represent an Rfam Hit to compute
    the QA information.
    """

    model = attr.ib(validator=is_a(basestring))
    model_rna_type = attr.ib(validator=is_a(basestring))
    model_domain = attr.ib(validator=optional(is_a(basestring)))
    model_completeness = attr.ib(validator=is_a(float), convert=float)
    sequence_completeness = attr.ib(validator=is_a(float), convert=float)

    @classmethod
    def build(cls, data):
        """
        Create a new RfamHit object. This accepts a dict where all keys match
        the attributes of this class.
        """
        return cls(**data)  # pylint: disable=star-args


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
    rfam_hits = attr.ib(validator=is_a(list))

    def is_species_specific(self):
        """
        Check if this sequence is specific to a species or is generic.
        """
        return self.taxid is not None

    def is_mitochondrial(self):
        """
        Check if this accession is mitochrondrial.
        """
        return any(a.is_mitochondrial() for a in self.accessions)

    def domains(self):
        """
        Get the set of all domains assigned to this sequence. This will ignore
        any invalid assignments (environmental samples, uncultured, or
        synthetic).
        """

        domains = set()
        for accession in self.accessions:
            domain = accession.domain
            if domain:
                domains.add(domain)
        return domains

    def has_unique_hit(self):
        """
        Check if there is only one Rfam hit to this sequence.
        """
        return len(self.rfam_hits) == 1


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
        hits = []
        for hit in data['hits']:
            if hit.pop('rfam_hit_id'):
                hits.append(RfamHit.build(hit))

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
            rfam_hits=hits,
        )


@attr.s()
class GenericSequence(Sequence):
    """
    This represents a sequence that is not assigned to a taxid.
    """

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
            rfam_hits=[],
        )
