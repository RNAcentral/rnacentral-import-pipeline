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

import hashlib

import pytest

from rfam import utils


def test_fetch_file_gets_unzipped_contents():
    contents = utils.fetch_file('12.2', 'database_files/family.txt.gz')
    assert hashlib.md5(contents.read()).hexdigest() == "5060f63e875739dab1ae514f5ae4b775"


def test_fetch_file_handles_uncompressed_files():
    contents = utils.fetch_file('12.2', 'USERMAN')
    assert hashlib.md5(contents.read()).hexdigest() == "6bb395c713e85a2db1b5d62ac96d4991"


@pytest.mark.parametrize("name,expected", [
    ("Cis52_sRNA", "other"),
    ("CbSR2", None),
    ("ppoRNA", "rRNA"),
    ("SNOR75", "snoRNA"),
    ("URE2_IRES", None),
    ("mir-393", "precursor_RNA"),
    ("GlsR20", "snoRNA"),
    ("Hammerhead_HH10", "hammerhead_ribozyme"),
    ("TTC28-AS1_3", "lncRNA"),
    ("ROSE_2", None),
    ("astro_FSE", None),
    ("RUF20", None),
    ("lactis-plasmid", "other"),
    ("leu-phe_leader", None),
    ("CRISPR-DR2", None),
    ('Metazoa_SRP', 'SRP_RNA'),
    ('TP53TG1_1', 'lncRNA'),
    ('group-II-D1D4-4', 'autocatalytically_spliced_intron'),
    ('Vault', 'vault_RNA'),
    ('tmRNA', 'tmRNA'),
    ('CDKN2B-AS', 'lncRNA'),
])
def test_can_fetch_a_mapping_from_name_to_isndc(name, expected):
    mapping = utils.name_to_insdc_type()
    assert mapping[name] == expected


@pytest.mark.parametrize("family_id,expected", [
    ("RF02647", "other"),
    ("RF02558", None),
    ("RF02554", "rRNA"),
])
def test_can_fetch_a_mapping_from_id_to_isndc(family_id, expected):
    mapping = utils.id_to_insdc_type()
    assert mapping[family_id] == expected
