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
import collections as coll

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from databases.data import Entry

CHROMOSOME_LEVEL_MAPPINGS = set([
    "WORMBASE",
    "POMBASE",
    "TAIR",
])

STANDARD_COLUMNS = set([
    'Source',
    'Target',
    'Source secondary accession',
    'Target secondary accession',
    'Source primary accession',
    'Target primary accession',
])


def tpa_key(value, database=None):
    """
    Generate a key that can be used to map from a Tpa to an Entry.
    """

    db_name = value.database
    if database:
        db_name = database

    if db_name in CHROMOSOME_LEVEL_MAPPINGS:
        if isinstance(value, Entry):
            if db_name == 'WORMBASE':
                return (value.parent_accession, value.standard_name)
            else:
                return (value.parent_accession, value.locus_tag)

        elif isinstance(value, Tpa):
            if db_name == 'POMBASE':
                return (value.parent_accession, value.database_accession)
            else:
                return (value.parent_accession, value.locus_tag)

    return (value.parent_accession, None)


def internal_database_name(ena_name):
    name = ena_name.upper()
    if name == 'POMBASE':
        return name
    if name == 'SGD':
        return name
    if name == 'SRPDB':
        return name
    if name == 'TAIR':
        return name
    if name == 'WORMBASE':
        return name
    if name == 'DICTYBASE':
        return name
    if name == 'LNCRNADB':
        return name
    if name == 'MIRBASE':
        return name
    if name == 'SNOPYDB':
        return 'SNOPY'
    if name == 'TMRNA-WEBSITE':
        return 'TMRNA_WEB'
    raise ValueError("Unknown database: %s" % ena_name)


@attr.s(frozen=True, slots=True)
class Tpa(object):
    database = attr.ib(validator=is_a(basestring))
    database_accession = attr.ib(validator=is_a(basestring))
    locus_tag = attr.ib(validator=optional(is_a(basestring)))
    parent_accession = attr.ib(validator=is_a(basestring))
    parent_secondary = attr.ib(validator=optional(is_a(basestring)))
    extra = attr.ib(
        validator=is_a(dict),
        default=attr.Factory(dict),
        hash=False,
    )

    @classmethod
    def from_dict(cls, row):
        locus_tag = row['Source secondary accession']
        if not locus_tag:
            locus_tag = None

        secondary = row['Target secondary accession']
        if not secondary:
            secondary = None

        extra = {k: v for k, v in row.items() if k not in STANDARD_COLUMNS}

        return cls(
            internal_database_name(row['Source']),
            row['Source primary accession'],
            locus_tag,
            row['Target primary accession'],
            secondary,
            extra,
        )


def tpa_attributes(tpa, entry):
    accession = '{accession}:{database}:{db_accession}'.format(
        accession=entry.accession,
        database=tpa.database,
        db_accession=tpa.database_accession,
    )

    secondary_id = None
    xrefs = entry.xref_data.get('ena_refs', {})
    if tpa.database in xrefs:
        _, secondary_id = xrefs[tpa.database]

    return {
        'primary_id': tpa.database_accession,
        'optional_id': secondary_id,
        'accession': accession,
        'database': tpa.database,
        'is_composite': 'Y',
        'non_coding_id': entry.accession,
    }


def silva_attributes(_, entry):
    if entry.database != 'SILVA':
        return {}
    return {}


def url_attributes(tpa, entry):
    """
    A class to build the url attributes on updated entries.
    """

    name = tpa.database.lower()
    if name == 'pombase':
        url = 'http://www.pombase.org/spombe/result/%s' % entry.primary_id

    elif name == 'sgd':
        url = 'http://www.yeastgenome.org/locus/%s/overview' % entry.primary_id

    elif name == 'srpdb':
        url = 'http://rnp.uthscsa.edu/rnp/SRPDB/rna/sequences/fasta/%s' % entry.primary_id

    elif name == 'tair':
        url = 'https://www.arabidopsis.org/servlets/TairObject?id=%s&type=locus' % entry.primary_id

    elif name == 'wormbase':
        url = 'http://www.wormbase.org/species/c_elegans/gene/%s' % entry.primary_id

    elif name == 'dictybase':
        url = 'http://dictybase.org/gene/%s' % entry.primary_id

    elif name == 'lncrnadb':
        url = 'http://www.lncrnadb.org/Detail.aspx?TKeyID=%s' % entry.primary_id

    elif name == 'mirbase':
        url = 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s' % entry.primary_id

    elif name == 'snopy':
        url = 'http://snoopy.med.miyazaki-u.ac.jp/snorna_db.cgi?mode=sno_info&id=%s' % entry.primary_id

    elif name == 'tmrna_web':
        url = 'http://bioinformatics.sandia.gov/tmrna/seqs/%s.html' % entry.primary_id

    else:
        raise ValueError("Unknown database for %s" % tpa)

    return {'url': url}


@attr.s(frozen=True)
class Transformer(object):
    attribute_builders = attr.ib()

    @classmethod
    def standard(cls):
        return cls([
            tpa_attributes,
            url_attributes,
            silva_attributes,
        ])

    def transform(self, tpa, entry):
        for builder in self.attribute_builders:
            entry = attr.evolve(entry, **builder(tpa, entry))
        return entry

    def transform_all(self, tpas, entry):
        if not tpas:
            return [entry]
        return [self.transform(t, entry) for t in tpas]


@attr.s()
class TpaMappings(object):
    databases = attr.ib(default=attr.Factory(set))
    simple_mapping = attr.ib(
        default=attr.Factory(lambda: coll.defaultdict(set))
    )

    def apply(self, entries):
        transformer = Transformer.standard()
        for entry in entries:
            tpas = self[entry]
            for transformed in transformer.transform_all(tpas, entry):
                yield transformed

    def update(self, tpas):
        for tpa in tpas:
            self.simple_mapping[tpa_key(tpa)].add(tpa)
            self.databases.add(tpa.database)

    def __getitem__(self, entry):
        for database in self.databases:
            key = tpa_key(entry, database=database)
            return self.simple_mapping.get(key, [])

    def __contains__(self, entry):
        return bool(self[entry])


def parse_tpa_file(handle):
    reader = csv.DictReader(handle, delimiter='\t')
    rows = it.ifilter(lambda r: r['Target'] == 'sequence', reader)
    return it.imap(Tpa.from_dict, rows)


def load(filenames):
    mapping = TpaMappings()
    for filename in filenames:
        with open(filename, 'rb') as raw:
            mapping.update(parse_tpa_file(raw))
    return mapping
