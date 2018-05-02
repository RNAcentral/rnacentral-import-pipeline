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
import itertools as it

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from ontologies import helpers as ont


def empty_to_none(raw):
    """
    Return an empty string to None otherwise return the string.
    """

    if not raw:
        return None
    return raw


@attr.s()
class RfamDatabaseLink(object):
    """
    This class represents the entries in the database_link in Rfam.
    """

    rfam_family = attr.ib(validator=is_a(basestring))
    database = attr.ib(validator=is_a(basestring))
    comment = attr.ib(
        convert=empty_to_none,
        validator=optional(is_a(basestring)),
    )
    external_id = attr.ib(validator=is_a(basestring))
    other = attr.ib(
        convert=empty_to_none,
        validator=optional(is_a(basestring)),
    )

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
        )

    def from_ontology(self):
        """
        Check if this instance comes from a known ontology (SO or GO).
        """
        return self.database in {'SO', 'GO'}

    def ontology_term(self):
        """
        Create an ontology term for this instance, if the instance comes from
        an ontology.
        """

        if not self.from_ontology():
            return None
        return ont.term(self.external_id)


def parse(handle):
    """
    Parse the given filehandle to produce all database link objects in the
    file.
    """

    reader = csv.DictReader(handle, delimiter='\t')
    return it.imap(RfamDatabaseLink.from_row, reader)


def ontology_terms(handle):
    """
    Produce an iterable of all ontology terms from Rfam.
    """

    seen = set()
    for reference in parse(handle):
        if not reference.from_ontology():
            continue
        if reference.external_id in seen:
            continue
        term = reference.ontology_term()
        seen.add(term.ontology_id)
        yield term
