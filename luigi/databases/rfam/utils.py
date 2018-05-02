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


def name_to_insdc_type(version='CURRENT'):
    families = load_families(version=version)
    return {family.name: family.guess_insdc() for family in families}


def id_to_insdc_type(version='CURRENT'):
    families = load_families(version=version)
    return {family.id: family.guess_insdc() for family in families}


def name_to_suppression(version='CURRENT'):
    families = load_families(version=version)
    return {family.name: family.is_suppressed for family in families}
