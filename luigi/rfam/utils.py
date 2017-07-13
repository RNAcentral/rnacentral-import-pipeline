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

import re
import csv
import gzip
from cStringIO import StringIO
from ftplib import FTP
import collections as coll

import attr
from attr.validators import instance_of as is_a
from functools32 import lru_cache

__doc__ = """
This module contains utility classes and functions for fetching and parsing
Rfam data from the FTP site.
"""


def as_0_based(raw):
    return int(raw) - 1


def convert_overlap(raw):
    if raw == '!' or raw == 'unique':
        return u'unique'
    if raw == '^' or raw == 'best':
        return u'best'
    raise Exception("Unknown overlap symbol %s" % raw)


def convert_strand(raw):
    if raw == '+' or raw == 1:
        return 1
    if raw == '-' or raw == -1:
        return -1
    raise Exception("Unknown strand %s" % raw)


def dash_none(raw):
    if raw == '-':
        return None
    return raw


def convert_trunc(raw):
    if isinstance(raw, tuple):
        return raw
    if raw == 'no':
        return (False, False)
    if raw == "5'":
        return (True, False)
    if raw == "3'":
        return (False, True)
    if raw == "5'&3'":
        return (True, True)
    raise Exception("Unknown yes/no %s" % raw)


def empty_str_from(target):
    def fn(raw):
        if raw == target:
            return u''
        return raw
    return fn


@attr.s(frozen=True, slots=True)
class RfamHit(object):
    target_name = attr.ib()
    seq_acc = attr.ib(convert=dash_none)
    rfam_name = attr.ib()
    rfam_acc = attr.ib()
    mdl = attr.ib()
    mdl_from = attr.ib(convert=as_0_based)
    mdl_to = attr.ib(convert=int)
    seq_from = attr.ib(convert=as_0_based)
    seq_to = attr.ib(convert=int)
    strand = attr.ib(convert=convert_strand)
    trunc = attr.ib(convert=convert_trunc)
    infernal_pass = attr.ib(convert=int)
    infernal_gc = attr.ib(convert=float)
    bias = attr.ib(convert=float)
    score = attr.ib(convert=float)
    e_value = attr.ib(convert=float)
    inc = attr.ib(convert=convert_overlap)
    description = attr.ib()


INFORMATIVE_NAMES = {
    # The leading _ is to prevent this from hitting isrP
    "_srp": "SRP_RNA",
    "y_rna": "Y_RNA",
    "hammerhead": "hammerhead_ribozyme",
    "group-ii": "autocatalytically_spliced_intron",
    "vault": "vault_RNA",
    "tmrna": "tmRNA",
    "rnase_p": "RNase_P_RNA",
    "rnasep": "RNase_P_RNA",
    "rnase_mrp": "RNase_MRP_RNA",
    "telomerase": "telomerase_RNA",

    # This is to prevent this from hitting 'ctRNA' as well
    "^trna": "tRNA",

    # These two patterns are so that we don't hit tracrRNA with rRNA
    "_rRNA": "rRNA",
    "PK-G12rRNA": "rRNA",
}

SO_TERM_MAPPING = {
    "SO:0000013": "scRNA",
    "SO:0000077": "antisense_RNA",
    "SO:0000252": "rRNA",
    "SO:0000253": "tRNA",
    "SO:0000274": "snRNA",
    "SO:0000275": "snoRNA",
    "SO:0000276": "precursor_RNA",
    "SO:0000370": "other",
    "SO:0000374": "ribozyme",
    "SO:0000380": "hammerhead_ribozyme",
    "SO:0000385": "RNase_MRP_RNA",
    "SO:0000386": "RNase_P_RNA",
    "SO:0000390": "telomerase_RNA",
    "SO:0000404": "vault_RNA",
    "SO:0000405": "Y_RNA",
    "SO:0000454": "rasiRNA",
    "SO:0000584": "tmRNA",
    "SO:0000588": "autocatalytically_spliced_intron",
    "SO:0000590": "SRP_RNA",
    "SO:0000602": "guide_RNA",
    "SO:0000644": "antisense_RNA",
    "SO:0000646": "siRNA",
    "SO:0000655": "other",
    "SO:0001035": "piRNA",
    "SO:0001244": "precursor_RNA",
    "SO:0001463": "lncRNA",
    "SO:0001641": "lncRNA",
    "SO:0001877": "lncRNA",
}

RFAM_RNA_TYPE_MAPPING = {
    "Cis-reg; leader;": "other",
    "Cis-reg; riboswitch;": "other",
    "Gene; CRISPR;": "other",
    "Gene; antisense;": "antisense_RNA",
    "Gene; antitoxin;": "other",
    "Gene; lncRNA;": "lncRNA",
    "Gene; miRNA;": "precursor_RNA",
    "Gene; rRNA;": "rRNA",
    "Gene; ribozyme;": "ribozyme",
    "Gene; sRNA;": "other",
    "Gene; snRNA; snoRNA; CD-box;": "snoRNA",
    "Gene; snRNA; snoRNA; HACA-box;": "snoRNA",
    "Gene; snRNA; snoRNA; scaRNA;": "snoRNA",
    "Gene; snRNA; splicing;": "snRNA",
    "Gene; snRNA;": "snRNA",
    "Gene; tRNA;": "tRNA",
    "Gene; tRNA;": "tRNA",
    "Gene;": "other",
}


@attr.s(frozen=True)
class RfamFamily(object):
    id = attr.ib(validator=is_a(basestring))
    name = attr.ib(validator=is_a(basestring))
    pretty_name = attr.ib(validator=is_a(basestring))
    so_terms = attr.ib(validator=is_a(set))
    rna_type = attr.ib(validator=is_a(basestring))
    domain = attr.ib()
    description = attr.ib(validator=is_a(basestring), convert=empty_str_from(r'\N'))
    seed_count = attr.ib(validator=is_a(int))
    full_count = attr.ib(validator=is_a(int))
    clan_id = attr.ib()
    length = attr.ib(validator=is_a(int))

    @classmethod
    def build_all(cls, clan_file, link_file, family_file):
        so_terms = coll.defaultdict(set)
        for line in link_file:
            parts = line.split()
            if parts[1] != 'SO':
                continue
            so_terms[parts[0]].add('SO:%s' % parts[2])

        clans = coll.defaultdict(set)
        for line in clan_file:
            parts = line.split()
            clans[parts[1]] = parts[0]

        families = []
        for row in csv.reader(family_file, delimiter='\t'):
            family = row[0]
            families.append(cls(
                id=family,
                name=row[1],
                pretty_name=row[3],
                so_terms=so_terms[family],
                rna_type=row[18].strip(),
                domain=None,
                description=row[9],
                seed_count=int(row[14]),
                full_count=int(row[15]),
                clan_id=clans.get(family, None),
                length=int(row[28]),
            ))
        return families

    @property
    def is_suppressed(self):
        return 'lncRNA' in self.rna_type or \
            not self.rna_type.startswith('Gene')

    def guess_insdc_using_name(self):
        found = set()
        for name, rna_type in INFORMATIVE_NAMES.items():
            if re.search(name, self.name, re.IGNORECASE):
                found.add(name)

        if found:
            if len(found) > 1:
                raise ValueError("Name patterns not distinct %s" %
                                 ', '.join(sorted(found)))
            return INFORMATIVE_NAMES[found.pop()]
        return None

    def guess_insdc_using_rna_type(self):
        return RFAM_RNA_TYPE_MAPPING.get(self.rna_type, None)

    def guess_insdc_using_so_terms(self):
        terms = set(SO_TERM_MAPPING.get(so, None) for so in self.so_terms)
        if len(terms) != 1:
            return None
        return terms.pop()

    def guess_insdc(self):
        return self.guess_insdc_using_name() or \
            self.guess_insdc_using_so_terms() or \
            self.guess_insdc_using_rna_type()


@attr.s()
class RfamClan(object):
    id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    description = attr.ib(validator=is_a(str))
    families = attr.ib(validator=is_a(set))

    @classmethod
    def build_all(cls, clan_file, membership_file):
        membership = coll.defaultdict(set)
        for line in membership_file:
            parts = line.split()
            membership[parts[0]].add(parts[1])
        membership = dict(membership)

        clans = []
        for line in clan_file:
            parts = line.strip().split('\t')
            clans.append(cls(
                id=parts[0],
                name=parts[3],
                description=parts[5],
                families=membership[parts[0]],
            ))
        return clans

    @property
    def family_count(self):
        return len(self.families)


def fetch_file(version, filename):
    ftp = FTP('ftp.ebi.ac.uk')
    ftp.login()
    ftp.cwd('pub/databases/Rfam/{version}'.format(
        version=version,
    ))
    command = 'RETR {filename}'.format(filename=filename)
    raw = StringIO()
    ftp.retrbinary(command, raw.write)
    ftp.quit()
    raw.seek(0)
    if filename.endswith('.gz'):
        raw = gzip.GzipFile(fileobj=raw, mode='rt')
    return raw


def get_family_file(version='CURRENT'):
    return fetch_file(version=version, filename='database_files/family.txt.gz')


def get_link_file(version='CURRENT'):
    return fetch_file(version=version,
                      filename='database_files/database_link.txt.gz')


def get_clan_membership(version='CURRENT'):
    return fetch_file(version=version,
                      filename='database_files/clan_membership.txt.gz')


def get_clans(version='CURRENT'):
    return fetch_file(version=version, filename='database_files/clan.txt.gz')


@lru_cache()
def load_families(version='CURRENT'):
    family_file = get_family_file(version=version)
    link_file = get_link_file(version=version)
    clan_file = get_clan_membership(version=version)
    return RfamFamily.build_all(clan_file, link_file, family_file)


@lru_cache()
def load_clans(version='CURRENT'):
    clan_file = get_clans(version=version)
    membership_file = get_clan_membership(version=version)
    return RfamClan.build_all(clan_file, membership_file)


def name_to_insdc_type(version='CURRENT'):
    families = load_families(version=version)
    return {family.name: family.guess_insdc() for family in families}


def id_to_insdc_type(version='CURRENT'):
    families = load_families(version=version)
    return {family.id: family.guess_insdc() for family in families}


def tbl_iterator(filename):
    with open(filename, 'rb') as raw:
        for line in raw:
            if line.startswith('#'):
                continue
            fields = attr.fields(RfamHit)
            parts = re.split(r'\s+', line.strip(), len(fields) - 1)
            yield RfamHit(*parts)


def name_to_suppression(version='CURRENT'):
    families = load_families(version=version)
    return {family.name: family.is_suppressed for family in families}
