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
            exons=[],
            rna_type='ncRNA',
            url='http://www.google.com',
            seq_version=None,
        )


def test_cannot_build_entry_with_empty_seq_id():
    with pytest.raises(TypeError):
        data.Entry(
            primary_id='a',
            accession='b',
            ncbi_tax_id=1,
            database='A',
            sequence='ACCG',
            exons=[],
            rna_type='ncRNA',
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
        exons=[],
        rna_type='ncRNA',
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
        exons=[],
        rna_type='ncRNA',
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.database == 'A_DATABASE_NAME'
    assert entry.database_name == 'A_DATABASE_NAME'


@pytest.mark.parametrize('rna_type', [
    'tRNA',
    'rRNA',
    'tmRNA',
    'misc_RNA',
    'precursor_RNA',
])
def test_labels_feature_type_rnas_correctly(rna_type):
    entry = data.Entry(
        primary_id='a',
        accession='b',
        ncbi_tax_id=1,
        database='a_database_name',
        sequence='ACCG',
        exons=[],
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
        exons=[],
        rna_type=rna_type,
        url='http://www.google.com',
        seq_version='1',
    )
    assert entry.feature_name == 'ncRNA'
    assert entry.ncrna_class is rna_type


def test_reference_can_handle_non_unicode():
    ref = data.Reference(
        authors=u'Ahn S., Jin T.E., Chang D.H., Rhee M.S., Kim H.J., Lee S.J., Park D.S., Kim B.C.',
        location='Int. J. Syst. Evol. Microbiol. 66(9):3656-3661(2016).',
        title=u'Agathobaculum butyriciproducens gen. nov. \xa0sp. nov., a strict anaerobic, butyrate-producing gut bacterium isolated from human faeces and reclassification of Eubacterium desmolans as Agathobaculum desmolans comb. nov',
        pmid=27334534,
        doi=None,
    )
    assert ref.md5() == '1d0712400fea453651afca2eb6ba0f4d'
