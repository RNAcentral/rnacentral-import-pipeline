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

import unittest as ut
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
    ('7SK', 'snRNA'),
    ('isrP', 'other'),
    ("U2", "snRNA"),
])
def test_can_fetch_a_mapping_from_name_to_isndc(name, expected):
    mapping = utils.name_to_insdc_type()
    assert mapping[name] == expected


@pytest.mark.parametrize("family_id,expected", [
    ("RF02647", "other"),
    ("RF02558", None),
    ("RF02554", "rRNA"),
    ("RF01268", "snoRNA"),
    ("RF00006", "vault_RNA"),
    ("RF01960", "rRNA"),
    ("RF01855", "SRP_RNA"),
    ("RF00017", "SRP_RNA"),
    ("RF02144", None),
    ("RF01414", None),
    ("RF00106", "antisense_RNA"),
    ("RF00242", "antisense_RNA"),
    ("RF01695", "antisense_RNA"),
    ("RF00039", "antisense_RNA"),
    ("RF00079", "other"),
    ("RF00003", "snRNA"),
    ("RF00548", "snRNA"),
    ("RF00024", "telomerase_RNA"),
    ("RF01050", "telomerase_RNA"),
    ("RF02277", "hammerhead_ribozyme"),
    ("RF01577", "RNase_P_RNA"),
    ("RF01807", "autocatalytically_spliced_intron"),
    ("RF02357", "RNase_P_RNA"),
    ("RF00094", "ribozyme"),
    ("RF02276", "hammerhead_ribozyme"),
    ("RF00023", "tmRNA"),
])
def test_can_fetch_a_mapping_from_id_to_isndc(family_id, expected):
    mapping = utils.id_to_insdc_type()
    assert mapping[family_id] == expected


def test_it_maps_a_known_lncRNA():
    mapping = utils.id_to_insdc_type()
    assert mapping['RF01800'] == 'lncRNA'


class INSDCRNATypeTest(ut.TestCase):

    def test_it_does_not_call_isrp_srp(self):
        rna = utils.RfamFamily(
            id='RF01398',
            name='isrP',
            so_terms=set(['SO:0001263']),
            rna_type='gene sRNA',
        )
        assert rna.guess_insdc_using_name() != 'SRP_RNA'
        assert rna.guess_insdc() != 'SRP_RNA'

    def test_it_does_not_call_ctrna_trna(self):
        rna = utils.RfamFamily(
            id='RF00236',
            name='ctRNA_pGA1',
            so_terms=set(['SO:0000644']),
            rna_type='Gene; antisense',
        )
        assert rna.guess_insdc_using_name() != 'tRNA'
        assert rna.guess_insdc() == 'antisense_RNA'

    def test_it_does_not_label_tracrrna_rrna(self):
        rna = utils.RfamFamily(
            id='RF02348',
            name='tracrRNA',
            so_terms=set(['SO:0000655']),
            rna_type='',
        )
        assert rna.guess_insdc_using_name() != 'tRNA'
        assert rna.guess_insdc() == 'other'
