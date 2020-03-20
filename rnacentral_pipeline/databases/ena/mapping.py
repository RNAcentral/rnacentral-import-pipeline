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
import collections as coll

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline.databases.data import Entry

CHROMOSOME_LEVEL_MAPPINGS = set([
    "WORMBASE",
])

DATABASES = {
    'SRPDB',
    'WormBase',
    'lncRNAdb',
    'snOPYdb',
    'tmRNA-Website',
}


def tpa_key(value, database=None):
    """
    Generate a key that can be used to map from a GenericTpa to an Entry.
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

        elif isinstance(value, GenericTpa):
            return (value.parent_accession, value.locus_tag)

    return (value.parent_accession, None)


def internal_database_name(ena_name):
    name = ena_name.upper()
    if name == 'SNOPYDB':
        return 'SNOPY'
    if name == 'TMRNA-WEBSITE':
        return 'TMRNA_WEB'
    return name


@attr.s(frozen=True, slots=True)
class GenericTpa(object):
    database = attr.ib(validator=is_a(str))
    database_accession = attr.ib(validator=is_a(str))
    locus_tag = attr.ib(validator=optional(is_a(str)))
    parent_accession = attr.ib(validator=is_a(str))
    parent_secondary = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def from_tsv(cls, row):
        locus_tag = row['Source secondary accession']
        if not locus_tag:
            locus_tag = None

        secondary = row['Target secondary accession']
        if not secondary:
            secondary = None

        database = internal_database_name(row['Source'])

        return cls(
            database,
            row['Source primary accession'],
            locus_tag,
            row['Target primary accession'],
            secondary,
        )

    def accession(self, entry):
        return '{accession}:{database}:{db_accession}'.format(
            accession=entry.accession,
            database=self.database,
            db_accession=self.database_accession,
        )

    def optional_id(self, entry):
        xrefs = entry.xref_data.get('ena_refs', {})
        if self.database in xrefs:
            return xrefs[self.database][1]
        return None

    def transform(self, entry):
        return attr.evolve(
            entry,
            primary_id=self.database_accession,
            optional_id=self.optional_id(entry),
            accession=self.accession(entry),
            database=self.database,
            is_composite='Y',
            non_coding_id=entry.accession,
        )


class UrlBuilder(object):
    def sgd(self, entry):
        return 'http://www.yeastgenome.org/locus/%s/overview' % entry.primary_id

    def srpdb(self, entry):
        return 'http://rnp.uthscsa.edu/rnp/SRPDB/rna/sequences/fasta/%s' % entry.primary_id

    def wormbase(self, entry):
        return 'http://www.wormbase.org/species/c_elegans/gene/%s' % entry.primary_id

    def dictybase(self, entry):
        return 'http://dictybase.org/gene/%s' % entry.primary_id

    def lncrnadb(self, entry):
        return 'http://www.lncrnadb.org/Detail.aspx?TKeyID=%s' % entry.primary_id

    def mirbase(self, entry):
        return 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s' % entry.primary_id

    def snopy(self, entry):
        return 'http://snoopy.med.miyazaki-u.ac.jp/snorna_db.cgi?mode=sno_info&id=%s' % entry.primary_id

    def tmrna_web(self, entry):
        return 'http://bioinformatics.sandia.gov/tmrna/seqs/%s.html' % entry.primary_id

    def transform(self, entry):
        builder = getattr(self, entry.database.lower())
        return attr.evolve(entry, url=builder(entry))


@attr.s()
class TpaMappings(object):
    databases: set = attr.ib(default=attr.Factory(set))
    simple_mapping = attr.ib(
        default=attr.Factory(lambda: coll.defaultdict(set))
    )
    counts = attr.ib(default=attr.Factory(coll.Counter))

    def add_tpas(self, tpas):
        for tpa in tpas:
            self.simple_mapping[tpa_key(tpa)].add(tpa)
            self.databases.add(tpa.database)
            self.counts[tpa.database] += 1

    def has_tpa_for(self, entry):
        return any(self.find_tpas(entry))

    def find_tpas(self, entry):
        for database in self.databases:
            key = tpa_key(entry, database=database)
            tpas = self.simple_mapping.get(key, [])
            for tpa in tpas:
                yield tpa
            if tpas:
                break

    def validate(self):
        dbs = [internal_database_name(db) for db in DATABASES]
        failed = [db for db in dbs if not self.counts[db]]
        if failed:
            raise ValueError("No TPAs found for: %s" % ', '.join(failed))
        return True


def parse_tpa_file(handle, klass=GenericTpa):
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        if row['Target'] != 'sequence':
            continue
        yield klass.from_tsv(row)


def load(raw):
    mapping = TpaMappings()
    mapping.add_tpas(parse_tpa_file(raw))
    return mapping


def load_file(filename):
    with open(filename, 'r') as raw:
        return load(raw)


def apply(mapping, entries):
    urls = UrlBuilder()
    for entry in entries:
        if mapping.has_tpa_for(entry):
            for tpa in mapping.find_tpas(entry):
                updated = tpa.transform(entry)
                yield urls.transform(updated)
        else:
            yield entry
