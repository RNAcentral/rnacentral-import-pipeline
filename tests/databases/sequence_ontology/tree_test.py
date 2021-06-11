# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import tempfile

from rnacentral_pipeline.databases.sequence_ontology import tree as so

import pytest


@pytest.fixture(scope='module')
def ontology():
    return so.load_ontology(so.REMOTE_ONTOLOGY)


@pytest.mark.parametrize(
    "so_term_id,expected",
    [
        ("rRNA", [("SO:0000655", "ncRNA"), ("SO:0000252", "rRNA")]),
        (
            "cytosolic_18S_rRNA",
            [
                ("SO:0000655", "ncRNA"),
                ("SO:0000252", "rRNA"),
                ('SO:0002343', 'cytosolic_rRNA'),
                ("SO:0000650", "cytosolic_SSU_rRNA"),
                ("SO:0000407", "cytosolic_18S_rRNA"),
            ],
        ),
        (
            "lincRNA",
            [
                ("SO:0000655", "ncRNA"),
                ("SO:0001877", "lncRNA"),
                ("SO:0001463", "lincRNA"),
            ],
        ),
        (
            "group_I_intron",
            [
                ("SO:0000188", "intron"),
                ("SO:0000588", "autocatalytically_spliced_intron"),
                ("SO:0000587", "group_I_intron"),
            ],
        ),
        (
            "ribozyme",
            [
                ("SO:0000673", "transcript"),
                ("SO:0000374", "ribozyme"),
            ],
        ),
        ("hammerhead_ribozyme", [
            ("SO:0000673", "transcript"),
            ("SO:0000374", "ribozyme"),
            ("SO:0000380", "hammerhead_ribozyme")
        ]),
        ("ncRNA", [("SO:0000655", "ncRNA")]),
        ("rRNA_primary_transcript", [
            ('SO:0000655', 'ncRNA'),
            ('SO:0000252', 'rRNA'),
            ('SO:0000655', 'rRNA_primary_transcript'),
        ]),
        ("antisense_lncRNA", [
            ("SO:0000655", "ncRNA"),
            ("SO:0001877", "lncRNA"),
            ("SO:0001904", "antisense_lncRNA"),
        ]),
        ("MicF_RNA", [
            ("SO:0000655", "ncRNA"),
            ("SO:0000644", "antisense_RNA"),
        ]),
        ("class_I_RNA", [
            ("SO:0000655", "ncRNA"),
        ]),
        ("RNA_6S", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("RprA_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("DsrA_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("CsrB_RsmB_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("spot_42_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("OxyS_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("RRE_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("GcvB_RNA", [
            ("SO:0000655", "ncRNA"),
            ('SO:0002247', 'sncRNA'),
            ("SO:0000370", "small_regulatory_ncRNA"),
        ]),
        ("pre_miRNA", [
            ("SO:0000673", "transcript"),
            ('SO:0001244', "pre_miRNA"),
        ]),
    ],
)
def test_can_compute_some_simple_paths(ontology, so_term_id, expected):
    assert so.rna_type_tree(ontology, so_term_id) == expected


@pytest.mark.parametrize("old_term,expected", [
    ("SO:0001171", [
        ("SO:0000655", "ncRNA"),
        ("SO:0000252", "rRNA"),
        ("SO:0002128", "mt_rRNA"),
        ("SO:0002345", "mt_LSU_rRNA"),
    ]),
    ("rRNA_21S", [
        ("SO:0000655", "ncRNA"),
        ("SO:0000252", "rRNA"),
        ("SO:0002128", "mt_rRNA"),
        ("SO:0002345", "mt_LSU_rRNA"),
    ]),
])
def test_can_compute_path_for_outdated_term(ontology, old_term, expected):
    tree =  so.rna_type_tree(ontology, old_term) 
    print(tree)
    assert tree == expected


@pytest.mark.parametrize("old_id,expected_name", [
    ("SO:0001171", "mt_LSU_rRNA"),
    ("SO:0002128", "mt_rRNA"),
])
def test_does_track_replacments(ontology, old_id, expected_name):
    assert ontology.id_to_name[old_id] == expected_name


@pytest.mark.parametrize(
    "so_id,name",
    [
        ("SO:0000584", "tmRNA"),
        ("SO:0000602", "guide_RNA"),
        ("SO:0000390", "telomerase_RNA"),
        ("SO:0001877", "lncRNA"),
    ],
)
def test_can_create_expected_mapping(ontology, so_id, name):
    with tempfile.NamedTemporaryFile() as tmp:
        mapping = so.name_index(ontology, tmp.name)
        assert mapping[so_id] == so_id
        assert mapping[name] == so_id
