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
import typing as ty

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

TpaKey = ty.Tuple[str, ty.Optional[str]]


def tpa_key(value: ty.Union[Entry, "GenericTpa"], database: ty.Optional[str]=None) -> TpaKey:
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


def internal_database_name(ena_name: str) -> str:
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
    def from_tsv(cls, row) -> "GenericTpa":
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

    def accession(self, entry: Entry) -> str:
        return f'{entry.accession}:{self.database}:{self.database_accession}'

    def optional_id(self, entry: Entry) -> ty.Optional[str]:
        xrefs = entry.xref_data.get('ena_refs', {})
        if self.database in xrefs:
            return xrefs[self.database][1]
        return None

    def transform(self, entry: Entry) -> Entry:
        return attr.evolve(
            entry,
            primary_id=self.database_accession,
            optional_id=self.optional_id(entry),
            accession=self.accession(entry),
            database=self.database,
            is_composite='Y',
            non_coding_id=entry.accession,
        )


class UrlBuilder:
    def sgd(self, entry) -> str:
        return 'http://www.yeastgenome.org/locus/%s/overview' % entry.primary_id

    def srpdb(self, entry: Entry) -> str:
        return 'http://rnp.uthscsa.edu/rnp/SRPDB/rna/sequences/fasta/%s' % entry.primary_id

    def wormbase(self, entry: Entry) -> str:
        return 'http://www.wormbase.org/species/c_elegans/gene/%s' % entry.primary_id

    def dictybase(self, entry: Entry) -> str:
        return 'http://dictybase.org/gene/%s' % entry.primary_id

    def lncrnadb(self, entry: Entry) -> str:
        return 'http://www.lncrnadb.org/Detail.aspx?TKeyID=%s' % entry.primary_id

    def mirbase(self, entry: Entry) -> str:
        return 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=%s' % entry.primary_id

    def snopy(self, entry: Entry) -> str:
        return 'http://snoopy.med.miyazaki-u.ac.jp/snorna_db.cgi?mode=sno_info&id=%s' % entry.primary_id

    def tmrna_web(self, entry: Entry) -> str:
        return 'http://bioinformatics.sandia.gov/tmrna/seqs/%s.html' % entry.primary_id

    def transform(self, entry: Entry) -> Entry:
        builder = getattr(self, entry.database.lower())
        return attr.evolve(entry, url=builder(entry))


@attr.s()
class TpaMappings:
    databases: set = attr.ib(default=attr.Factory(set))
    simple_mapping: ty.Dict[TpaKey, ty.Set[GenericTpa]] = attr.ib(
        default=attr.Factory(lambda: coll.defaultdict(set))
    )
    counts: ty.Dict[str, int] = attr.ib(default=attr.Factory(coll.Counter))

    def add_tpas(self, tpas: ty.Iterable[GenericTpa]):
        for tpa in tpas:
            self.simple_mapping[tpa_key(tpa)].add(tpa)
            self.databases.add(tpa.database)
            self.counts[tpa.database] += 1

    def has_tpa_for(self, entry: Entry) -> bool:
        return any(self.find_tpas(entry))

    def find_tpas(self, entry: Entry) -> ty.Iterable[GenericTpa]:
        for database in self.databases:
            key = tpa_key(entry, database=database)
            tpas = self.simple_mapping.get(key, set())
            for tpa in tpas:
                yield tpa
            if tpas:
                break

    def validate(self) -> bool:
        dbs = [internal_database_name(db) for db in DATABASES]
        failed = [db for db in dbs if not self.counts[db]]
        if failed:
            raise ValueError("No TPAs found for: %s" % ', '.join(failed))
        return True


def parse_tpa_file(handle) -> ty.Iterable[GenericTpa]:
    reader = csv.DictReader(handle, delimiter='\t')
    for row in reader:
        if row['Target'] != 'sequence':
            continue
        yield GenericTpa.from_tsv(row)


def load(raw) -> TpaMappings:
    mapping = TpaMappings()
    mapping.add_tpas(parse_tpa_file(raw))
    return mapping


def load_file(filename: str) -> TpaMappings:
    with open(filename, 'r') as raw:
        return load(raw)


def apply(mapping: TpaMappings, entries: ty.Iterable[Entry]) -> ty.Iterable[Entry]:
    urls = UrlBuilder()
    for entry in entries:
        if mapping.has_tpa_for(entry):
            for tpa in mapping.find_tpas(entry):
                updated = tpa.transform(entry)
                yield urls.transform(updated)
        else:
            yield entry
