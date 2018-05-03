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


EXCLUDED_TERMS = {
    'GO:0008049',
    'GO:0042981',
    'GO:0042749',
    'GO:0050789',
    'GO:0006810',
}

EXCLUDED_MIRNA = {'GO:0035068', 'GO:0006396'}


GO_REPLACEMENTS = {
    'GO:0044005': 'GO:0051819',
}


QUERY = """
select distinct
	link.*,
	family.type
from database_link link
join family on family.rfam_acc = link.rfam_acc
;
"""


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
    family_type = attr.ib(validator=is_a(basestring))

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


def correct_go_term(reference):
    """
    This will correct the reference to GO if required. Basically, this will
    exclude some terms and replace others.
    """

    go_term_id = reference.external_id
    if go_term_id in EXCLUDED_TERMS:
        return None

    if reference.rna_type == 'Gene; miRNA;' and go_term_id in EXCLUDED_MIRNA:
        return None

    go_term_id = GO_REPLACEMENTS.get(go_term_id, go_term_id)
    return attr.assoc(reference, external_id=go_term_id)


def ontology_references(handle):
    """
    Produce an iterable of all ontology terms from Rfam.
    """

    for reference in parse(handle):
        if not reference.from_ontology():
            continue

        if reference.database == 'GO':
            reference = correct_go_term(reference)

        if not reference:
            continue

        yield reference


def ontology_terms(handle):
    seen = set()
    for reference in ontology_references(handle):
        if reference.external_id in seen:
            continue

        term = reference.ontology_term()
        seen.add(term.ontology_id)
        yield term
