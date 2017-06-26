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
import collections as coll

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
    ("CbSR2", "antisense_RNA"),
    ("ppoRNA", "rRNA"),
    ("SNOR75", "snoRNA"),
    ("URE2_IRES", None),
    ("mir-393", "precursor_RNA"),
    ("GlsR20", "snoRNA"),
    ("Hammerhead_HH10", "hammerhead_ribozyme"),
    ("TTC28-AS1_3", "lncRNA"),
    ("ROSE_2", None),
    ("astro_FSE", None),
    ("RUF20", 'other'),
    ("lactis-plasmid", "other"),
    ("leu-phe_leader", 'other'),
    ("CRISPR-DR2", 'other'),
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
    ("RF00003", "snRNA"),
    ("RF00006", "vault_RNA"),
    ("RF00008", "hammerhead_ribozyme"),
    ("RF00009", "RNase_P_RNA"),
    ("RF00010", "RNase_P_RNA"),
    ("RF00011", "RNase_P_RNA"),
    ("RF00017", "SRP_RNA"),
    ("RF00017", "SRP_RNA"),
    ("RF00023", "tmRNA"),
    ("RF00024", "telomerase_RNA"),
    ("RF00025", "telomerase_RNA"),
    ("RF00030", "RNase_MRP_RNA"),
    ("RF00039", "antisense_RNA"),
    ("RF00079", "other"),
    ("RF00094", "ribozyme"),
    ("RF00106", "antisense_RNA"),
    ("RF00163", "hammerhead_ribozyme"),
    ("RF00169", "SRP_RNA"),
    ("RF00242", "antisense_RNA"),
    ("RF00373", "RNase_P_RNA"),
    ("RF00548", "snRNA"),
    ("RF01050", "telomerase_RNA"),
    ("RF01268", "snoRNA"),
    ("RF01414", 'other'),
    ("RF01502", "SRP_RNA"),
    ("RF01577", "RNase_P_RNA"),
    ("RF01695", "antisense_RNA"),
    ("RF01807", "autocatalytically_spliced_intron"),
    ("RF01854", "SRP_RNA"),
    ("RF01855", "SRP_RNA"),
    ("RF01856", "SRP_RNA"),
    ("RF01857", "SRP_RNA"),
    ("RF01960", "rRNA"),
    ("RF02144", 'other'),
    ("RF02276", "hammerhead_ribozyme"),
    ("RF02277", "hammerhead_ribozyme"),
    ("RF02357", "RNase_P_RNA"),
    ("RF02554", "rRNA"),
    ("RF02558", "antisense_RNA"),
    ("RF02647", "other"),
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
            pretty_name='Hfq binding RNA',
            so_terms=set(['SO:0001263']),
            rna_type='gene sRNA',
            domain=None,
            description='',
            seed_count=-1,
            full_count=-1,
            clan_id=None,
            length=148,
        )
        assert rna.guess_insdc_using_name() != 'SRP_RNA'
        assert rna.guess_insdc() != 'SRP_RNA'

    def test_it_does_not_call_ctrna_trna(self):
        rna = utils.RfamFamily(
            id='RF00236',
            name='ctRNA_pGA1',
            pretty_name='ctRNA',
            so_terms=set(['SO:0000644']),
            rna_type='Gene; antisense',
            domain=None,
            description='',
            seed_count=-1,
            full_count=-1,
            clan_id=None,
            length=79,
        )
        assert rna.guess_insdc_using_name() != 'tRNA'
        assert rna.guess_insdc() == 'antisense_RNA'

    def test_it_does_not_label_tracrrna_rrna(self):
        rna = utils.RfamFamily(
            id='RF02348',
            name='tracrRNA',
            pretty_name='Trans-activating crRNA',
            so_terms=set(['SO:0000655']),
            rna_type='',
            domain=None,
            description='',
            seed_count=-1,
            full_count=-1,
            clan_id=None,
            length=91,
        )
        assert rna.guess_insdc_using_name() != 'tRNA'
        assert rna.guess_insdc() == 'other'


class LoadingFamiliesTest(ut.TestCase):
    def test_it_loads_all_families(self):
        assert len(utils.load_families()) == 2588

    def test_it_can_load_family_correctly(self):
        assert utils.load_families()[0] == utils.RfamFamily(
            id='RF00001',
            name='5S_rRNA',
            pretty_name='5S ribosomal RNA',
            so_terms=set(['SO:0000652']),
            rna_type='Gene; rRNA;',
            domain=None,
            description=(
                '5S ribosomal RNA (5S rRNA) is a component of the '
                'large ribosomal subunit in both prokaryotes and eukaryotes. '
                'In eukaryotes, it is synthesised by RNA polymerase III (the '
                'other eukaryotic rRNAs are cleaved from a 45S precursor '
                'synthesised by RNA polymerase I). In Xenopus oocytes, it has '
                'been shown that fingers 4-7 of the nine-zinc finger '
                'transcription factor TFIIIA can bind to the central region '
                'of 5S RNA. Thus, in addition to positively regulating 5S '
                'rRNA transcription, TFIIIA also stabilises 5S rRNA until it '
                'is required for transcription.'
            ),
            seed_count=712,
            full_count=183439,
            clan_id=None,
            length=119,
        )

    def test_it_can_assign_correct_clan_ids(self):
        clans = coll.defaultdict(set)
        for f in utils.load_families():
            clans[f.clan_id].add(f.id)

        assert len(clans) == 111
        assert clans['CL00001'] == set([
            'RF00005',
            'RF00023',
            'RF01849',
            'RF01850',
            'RF01851',
            'RF01852',
            'RF02544',
        ])

    @pytest.mark.skip()
    def test_it_can_correctly_find_all_bacterial_families(self):
        families = utils.load_families()
        bacterial = set(f.id for f in families if f.domain == 'Bacteria')
        assert bacterial == set([
        ])

    def test_it_can_set_description_to_empty(self):
        families = utils.load_families()
        family = next(f for f in families if f.id == 'RF02493')
        assert family.description == ''


class LoadingClansTest(ut.TestCase):
    def test_it_can_load_all_clans(self):
        assert len(utils.load_clans()) == 110

    def test_it_can_load_a_clan_correctly(self):
        assert utils.load_clans()[0] == utils.RfamClan(
            id='CL00001',
            name='tRNA clan',
            description=(
                'The tRNA clan contains the RNA families tRNA and '
                'tmRNA. Homology between these families has been established '
                'in the published literature [1-5].'
            ),
            families=set([
                'RF00005',
                'RF00023',
                'RF01849',
                'RF01850',
                'RF01851',
                'RF01852',
                'RF02544',
            ])
        )
