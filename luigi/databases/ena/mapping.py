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
import operator as op

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

CHROMOSOME_LEVEL_MAPPINGS = set([
    "WORMBASE",
])


def tpa_key(value):
    if value.database in CHROMOSOME_LEVEL_MAPPINGS:
        return (value.parent_accession, value.locus_tag)
    return (value.parent_accession, None)


@attr.s(frozen=True, slots=True)
class GenericTpa(object):
    database = attr.ib(validator=is_a(basestring))
    database_accession = attr.ib(validator=is_a(basestring))
    locus_tag = attr.ib(validator=optional(is_a(basestring)))
    parent_accession = attr.ib(validator=is_a(basestring))

    @classmethod
    def from_tsv(cls, row):
        secondary = row['Source secondary accession']
        if not secondary:
            secondary = None

        database = row['Source'].upper()
        if database == 'SNOPYDB':
            database = 'SNOPY'
        if database == 'TMRNA-WEBSITE':
            database = 'TMRNA_WEB'

        return cls(
            database,
            row['Source primary accession'],
            secondary,
            row['Target primary accession'],
        )

    def accession(self, entry):
        return '{accession}:{database}:{db_accession}'.format(
            accession=entry.accession,
            database=self.database,
            db_accession=self.database_accession
        )

    def transform(self, entry):
        return attr.assoc(
            entry,
            primary_id=self.database_accession,
            accession=self.accession(entry),
            database=self.database,
            is_composite='Y',
            non_coding_id=entry.accession,
        )


class UrlBuilder(object):
    def pombase(self, entry):
        return 'http://www.pombase.org/spombe/result/%s' % entry.primary_id

    def sgd(self, entry):
        return 'http://www.yeastgenome.org/locus/%s/overview' % entry.primary_id

    def srpdb(self, entry):
        return 'http://rnp.uthscsa.edu/rnp/SRPDB/rna/sequences/fasta/%s' % entry.primary_id

    def tair(self, entry):
        return 'https://www.arabidopsis.org/servlets/TairObject?id=%s&type=locus' % entry.primary_id

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
        return attr.assoc(entry, url=builder(entry))


@attr.s()
class TpaMappings(object):
    simple_mapping = attr.ib(
        default=attr.Factory(lambda: coll.defaultdict(set))
    )

    def add_tpas(self, tpas):
        for tpa in tpas:
            self.simple_mapping[tpa_key(tpa)].add(tpa)

    def has_tpa_for(self, entry):
        return any(self.find_tpas(entry))

    def find_tpas(self, entry):
        for tpa in self.simple_mapping.get(tpa_key(entry), []):
            yield tpa


def parse_tpa_file(handle, klass=GenericTpa):
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        if row['Target'] != 'sequence':
            continue
        yield klass.from_tsv(row)


def load(filenames):
    mapping = TpaMappings()
    for filename in filenames:
        with open(filename, 'rb') as raw:
            mapping.add_tpas(parse_tpa_file(raw))
    return mapping


def apply(mapping, entries):
    urls = UrlBuilder()
    for entry in entries:
        if mapping.has_tpa_for(entry):
            for tpa in mapping.find_tpas(entry):
                updated = tpa.transform(entry)
                yield urls.transform(updated)
        else:
            yield entry
