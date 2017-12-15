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


@attr.s(frozen=True)
class GenericTpa(object):
    database = attr.ib(validator=is_a(basestring))
    database_accession = attr.ib(validator=is_a(basestring))
    locus_tag = attr.ib(validator=is_a(basestring))
    parent_accession = attr.ib(validator=is_a(basestring))
    feature_type = attr.ib(validator=is_a(basestring))

    @classmethod
    def from_tsv(cls, row):
        return cls(
            row['Source'].upper(),
            row['Source primary accession'],
            row['Source secondary accession'],
            row['Target secondary accession'],
            row['Target'],
        )

    def accession(self, entry):
        return '{accession}:{database}:{db_accession}'.format(
            entry.accession,
            self.database,
            self.database_accession
        )

    def transform(self, entry):
        return attr.assoc(
            entry,
            accession=self.accession(entry),
            database=self.database,
            is_composite='Y',
            non_coding_id=entry.accession,
        )


class UrlBuilder(object):
    def wormbase(self, entry):
        return 'http://www.wormbase.org/species/c_elegans/gene/%s' % entry.primary_id

    def srpdb(self, entry):
        return 'http://rnp.uthscsa.edu/rnp/SRPDB/rna/sequences/fasta/%s' % entry.primary_id

    def mirbase(self, entry):
        return 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s' % entry.primary_id

    def tmrna_website(self, entry):
        return 'http://bioinformatics.sandia.gov/tmrna/seqs/%s.html' % entry.primary_id

    def lncrnadb(self, entry):
        return 'http://www.lncrnadb.org/Detail.aspx?TKeyID=%s' % entry.primary_id

    def snopydb(self, entry):
        return 'http://snoopy.med.miyazaki-u.ac.jp/snorna_db.cgi?mode=sno_info&id=%s' % entry.primary_id

    def transform(self, entry):
        url_builder = getattr(self, entry.database.lower())
        return attr.assoc(entry, url=url_builder())


@attr.s()
class TpaMappings(object):
    simple_mapping = attr.ib(
        validator=is_a(dict),
        default=lambda: coll.defaultdict(set)
    )

    def add_tpas(self, tpas):
        for tpa in tpas:
            self.simple_mapping[tpas.parent_accession].add(tpa)

    def has_tpa_for(self, entry):
        return any(self.find_tpas(entry))

    def find_tpas(self, entry):
        for tpa in self.simple_mapping.get(entry.parent_accession, []):
            yield tpa


def parse_tpa_file(filename, klass):
    with open(filename, 'rb') as raw:
        reader = csv.DictReader(raw, delimiter='\t')
        for row in reader:
            if row['Target'] != 'sequence':
                continue
            yield klass.from_tsv(row)


def load(description):
    mapping = TpaMappings()
    for _, filename in description:
        with open(filename, 'rb') as raw:
            mapping.add_tpas(parse_tpa_file(raw, GenericTpa))
    return mapping


def apply(mapping, entries):
    urls = UrlBuilder()
    for entry in entries:
        if mapping.has_tpa_for(entry):
            for tpa in mapping.find_tpas(entry):
                updated = tpa.transform(entry)
                updated = urls.transform(updated)
                yield updated
        else:
            yield entry


# @attr.s()
# class WormBaseTpa(object):
#     wormbase_accession = attr.ib(validator=is_a(basestring))
#     locus_tag = attr.ib(validator=is_a(basestring))
#     ena_accession = attr.ib(validator=is_a(basestring))
#     feature_type = attr.ib(validator=is_a(basestring))
#     @classmethod
#     def from_tsv(cls, row):
#         return cls(
#             row['Source primary accession'],
#             row['Source secondary accession'],
#             row['Target secondary accession'],
#             row['Target'],
#         )
#     def accession(self, entry):
#         return '%s:WORMBASE:%s' % (entry.accession, self.wormbase_accession)
#     def transform(self, entry, extra={}):
#         updated = attr.assoc(
#             entry,
#             accession=self.accession(entry),
#             database='WORMBASE',
#             is_composite='Y',
#             non_coding_id=entry.accession,
#         )
#         return attr.assoc(updated, **extra)
