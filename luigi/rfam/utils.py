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

__doc__ = """
This module contains utility classes and functions for fetching and parsing
Rfam data from the FTP site.
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

INFORMATIVE_NAMES = {
    "srp": "SRP_RNA",
    "y_rna": "Y_RNA",
    "hammerhead": "hammerhead_ribozyme",
    "group-ii": "autocatalytically_spliced_intron",
    "vault": "vault_RNA",
    "tmrna": "tmRNA",
    "rnase_p": "RNase_P_RNA",
    "rnase_mrp": "RNase_MRP_RNA",
    "telomerase": "telomerase_RNA"
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
    "Gene; lncRNA;": "lncRNA",
    "Gene; miRNA;": "precursor_RNA",
    "Gene; rRNA;": "rRNA",
    "Gene; sRNA;": "other",
    "Gene; snRNA;": "snRNA",
    "Gene; snRNA; snoRNA; CD-box;": "snoRNA",
    "Gene; snRNA; snoRNA; HACA-box;": "snoRNA",
    "Gene; snRNA; snoRNA; scaRNA;": "snoRNA",
    "Gene; tRNA;": "tRNA"
}


@attr.s(frozen=True)
class RfamFamily(object):
    id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    so_terms = attr.ib(validator=is_a(set))
    rna_type = attr.ib(validator=is_a(str))

    @classmethod
    def build_all(cls, link_file, family_file):
        so_terms = coll.defaultdict(set)
        for line in link_file:
            parts = line.split()
            if parts[1] != 'SO':
                continue
            so_terms[parts[0]].add('SO:%s' % parts[2])

        families = []
        for row in csv.reader(family_file, delimiter='\t'):
            family = row[0]
            families.append(cls(
                id=family,
                name=row[1],
                so_terms=so_terms[family],
                rna_type=row[18].strip()
            ))
        return families

    def guess_insdc_using_name(self):
        for name, rna_type in INFORMATIVE_NAMES.items():
            if re.search(name, self.name, re.IGNORECASE):
                return rna_type
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
        raw = gzip.GzipFile(fileobj=raw, mode='rb')
    return raw


def get_families(version='CURRENT'):
    return fetch_file(version=version, filename='database_files/family.txt.gz')


def get_links(version='CURRENT'):
    return fetch_file(version=version,
                      filename='database_files/database_link.txt.gz')


@lru_cache()
def name_to_insdc_type(version='CURRENT'):
    family_file = get_families(version=version)
    link_file = get_links(version=version)
    families = RfamFamily.build_all(link_file, family_file)
    return {family.name: family.guess_insdc() for family in families}


@lru_cache()
def id_to_insdc_type(version='CURRENT'):
    family_file = get_families(version=version)
    link_file = get_links(version=version)
    families = RfamFamily.build_all(link_file, family_file)
    return {family.id: family.guess_insdc() for family in families}
