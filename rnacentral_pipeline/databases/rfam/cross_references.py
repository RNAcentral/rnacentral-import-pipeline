# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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
import operator as op
import itertools as it

import six

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.writers import MultiCsvOutput


EXCLUDED_TERMS = {
    'GO:0008049',
    'GO:0042981',
    'GO:0042749',
    'GO:0050789',
    'GO:0006810',
    'GO:0001263',
    'SO:0004725',
    'SO:0010039',
}

EXCLUDED_MIRNA = {'GO:0035068', 'GO:0006396'}


GO_REPLACEMENTS = {
    'GO:0044005': 'GO:0051819',
}


def empty_to_none(raw):
    """
    Return an empty string to None otherwise return the string.
    """

    if not raw:
        return None
    return raw


@attr.s(frozen=True)
class RfamDatabaseLink(object):
    """
    This class represents the entries in the database_link in Rfam.
    """

    rfam_family = attr.ib(validator=is_a(str))
    database = attr.ib(validator=is_a(str))
    comment = attr.ib(
        converter=empty_to_none,
        validator=optional(is_a(str)),
    )
    external_id = attr.ib(validator=is_a(str))
    other = attr.ib(
        converter=empty_to_none,
        validator=optional(is_a(str)),
    )
    family_type = attr.ib(validator=is_a(str))

    @classmethod
    def from_row(cls, row):
        """
        Build an object from a dictionary that has the same names as the
        database columns.
        """

        database = row['db_id']
        external_id = row['db_link']
        if database in {'SO', 'GO'}:
            external_id = '%s:%s' % (database, row['db_link'])

        return cls(
            rfam_family=row['rfam_acc'],
            database=database,
            comment=row['comment'],
            external_id=external_id,
            other=row['other_params'],
            family_type=row['type'],
        )

    def from_ontology(self):
        """
        Check if this instance comes from a known ontology (SO or GO).
        """
        return self.database in {'SO', 'GO'}

    def writeable_go_mappings(self):
        if self.database != 'GO':
            return

        yield [
            self.rfam_family,
            self.external_id,
        ]

    def writeable_ontology_terms(self):
        if not self.from_ontology():
            return

        yield [self.external_id]


def parse(handle):
    """
    Parse the given filehandle to produce all database link objects in the
    file.
    """

    reader = csv.DictReader(handle, delimiter='\t')
    return six.moves.map(RfamDatabaseLink.from_row, reader)


def correct_go_term(reference):
    """
    This will correct the reference to GO if required. Basically, this will
    exclude some terms and replace others.
    """

    go_term_id = reference.external_id
    if go_term_id in EXCLUDED_TERMS:
        return None

    if reference.family_type == 'Gene; miRNA;' and \
            go_term_id in EXCLUDED_MIRNA:
        return None

    go_term_id = GO_REPLACEMENTS.get(go_term_id, go_term_id)
    return attr.evolve(reference, external_id=go_term_id)


def ontology_references(handle):
    """
    Produce an iterable of all ontology terms from Rfam.
    """

    for reference in parse(handle):
        if not reference.from_ontology():
            continue

        if reference.database in {'SO', 'GO'}:
            reference = correct_go_term(reference)

        if not reference:
            continue

        yield reference


def from_file(handle, output):
    writer = MultiCsvOutput.build(
        ontology_references,
        terms={
            'transformer': op.methodcaller('writeable_ontology_terms'),
        },
        rfam_ontology_mappings={
            'transformer': op.methodcaller('writeable_go_mappings'),
        }
    )
    writer(output, handle)
