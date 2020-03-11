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

import pytest

from rnacentral_pipeline.databases import data


def test_cannot_build_entry_without_seq_id():
    with pytest.raises(TypeError):
        data.Entry(
            primary_id='a',
            accession='b',
            ncbi_tax_id=1,
            database='A',
            sequence='ACCG',
            regions=[],
            rna_type='SO:0000655',
            url='http://www.google.com',
            seq_version=None,
        )


def test_cannot_build_entry_with_empty_seq_id():
    with pytest.raises(ValueError):
        data.Entry(
            primary_id='a',
            accession='b',
            ncbi_tax_id=1,
            database='A',
            sequence='ACCG',
            regions=[],
            rna_type='SO:0000655',
            url='http://www.google.com',
            seq_version='',
        )


def test_can_build_with_seq_version():
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='A',
        sequence='ACCG',
        regions=[],
        rna_type='SO:0000655',
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.seq_version == '1'


def test_it_will_always_uppercase_database():
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='a_database_name',
        sequence='ACCG',
        regions=[],
        rna_type='SO:0000655',
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.database == 'A_DATABASE_NAME'
    assert entry.database_name == 'A_DATABASE_NAME'


@pytest.mark.parametrize('so_term,rna_type', [
    ('SO:0000385', 'RNase_MRP_RNA'),
    ('SO:0000386', 'RNase_P_RNA'),
    ('SO:0000590', 'SRP_RNA'),
    ('SO:0000405', 'Y_RNA'),
    ('SO:0000644', 'antisense_RNA'),
    ('SO:0000588', 'autocatalytically_spliced_intron'),
    ('SO:0000602', 'guide_RNA'),
    ('SO:0000380', 'hammerhead_ribozyme'),
    ('SO:0001877', 'lncRNA'),
    ('SO:0000276', 'miRNA'),
    ('SO:0000655', 'other'),
    ('SO:0001035', 'piRNA'),
    ('SO:0000454', 'rasiRNA'),
    ('SO:0000374', 'ribozyme'),
    ('SO:0000013', 'scRNA'),
    ('SO:0000646', 'siRNA'),
    ('SO:0000274', 'snRNA'),
    ('SO:0000275', 'snoRNA'),
    ('SO:0000390', 'telomerase_RNA'),
    ('SO:0000404', 'vault_RNA'),
    ('SO:0000253', 'tRNA'),
    ('SO:0000252', 'rRNA'),
    ('SO:0000584', 'tmRNA'),
    ('SO:0000673', 'misc_RNA'),
    # 'precursor_RNA',
])
def test_convertes_rna_types_correctl(rna_type, so_term):
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='a_database_name',
        sequence='ACCG',
        regions=[],
        rna_type=rna_type,
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.rna_type == so_term


@pytest.mark.parametrize('rna_type', [
    'tRNA',
    'rRNA',
    'tmRNA',
    'misc_RNA',
    # 'precursor_RNA',
])
def test_labels_feature_type_rnas_correctly(rna_type):
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='a_database_name',
        sequence='ACCG',
        regions=[],
        rna_type=rna_type,
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.feature_name == rna_type
    assert entry.ncrna_class is None


@pytest.mark.parametrize('rna_type', [
    "RNase_MRP_RNA",
    "RNase_P_RNA",
    "SRP_RNA",
    "Y_RNA",
    "antisense_RNA",
    "autocatalytically_spliced_intron",
    "guide_RNA",
    "hammerhead_ribozyme",
    "lncRNA",
    "miRNA",
    "other",
    "piRNA",
    "rasiRNA",
    "ribozyme",
    "scRNA",
    "siRNA",
    "snRNA",
    "snoRNA",
    "telomerase_RNA",
    "vault_RNA",
])
def test_labels_ncrna_types_correctly(rna_type):
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='a_database_name',
        sequence='ACCG',
        regions=[],
        rna_type=rna_type,
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.feature_name == 'ncRNA'
    assert entry.ncrna_class is rna_type


def test_can_write_valid_sequence_regions():
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='a_database_name',
        sequence='ACCGGGGGGGGGGGGGGGGGGGGGGGG',
        regions=[
            data.SequenceRegion(
                chromosome='1',
                strand=-1,
                exons=[
                    data.Exon(start=1, stop=10),
                    data.Exon(start=12, stop=14),
                ],
                assembly_id='hg19',
                coordinate_system=data.CoordinateSystem.one_based(),
            ),
            data.SequenceRegion(
                chromosome='12',
                strand=1,
                exons=[
                    data.Exon(start=20, stop=22),
                    data.Exon(start=24, stop=33),
                ],
                assembly_id='hg38',
                coordinate_system=data.CoordinateSystem.one_based(),
            )
        ],
        rna_type='snoRNA',
        url='http://www.google.com',
        seq_version='1',
    )

    assert list(entry.write_sequence_regions()) == [
        ['b', '@1/1-10,12-14:-', '1', -1, 'hg19', 2, 1, 10],
        ['b', '@1/1-10,12-14:-', '1', -1, 'hg19', 2, 12, 14],
        ['b', '@12/20-22,24-33:+', '12', 1, 'hg38', 2, 20, 22],
        ['b', '@12/20-22,24-33:+', '12', 1, 'hg38', 2, 24, 33],
    ]
