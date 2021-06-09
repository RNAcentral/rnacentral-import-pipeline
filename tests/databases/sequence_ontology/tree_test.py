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
from functools import lru_cache

from rnacentral_pipeline.databases.sequence_ontology import tree as so

import pytest


@lru_cache()
def load():
    return so.load_ontology(so.REMOTE_ONTOLOGY)


@pytest.mark.parametrize(
    "so_term_id,expected",
    [
        ("rRNA", [("SO:0000655", "ncRNA"), ("SO:0000252", "rRNA")]),
        (
            "rRNA_18S",
            [
                ("SO:0000655", "ncRNA"),
                ("SO:0000252", "rRNA"),
                ("SO:0000650", "small_subunit_rRNA"),
                ("SO:0000407", "rRNA_18S"),
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
                ("SO:0000372", "enzymatic_RNA"),
                ("SO:0000374", "ribozyme"),
            ],
        ),
        ("hammerhead_ribozyme", [("SO:0000374", "ribozyme"), ("SO:0000380", "hammerhead_ribozyme")]),
        ("ncRNA", [("SO:0000655", "ncRNA")]),
    ],
)
def test_can_compute_some_simple_paths(so_term_id, expected):
    data = load()
    assert so.rna_type_tree(data, so_term_id) == expected


@pytest.mark.parametrize(
    "so_id,name",
    [
        ("SO:0000584", "tmRNA"),
        ("SO:0000602", "guide_RNA"),
        ("SO:0000390", "telomerase_RNA"),
        ("SO:0001877", "lncRNA"),
    ],
)
def test_can_create_expected_mapping(so_id, name):
    with tempfile.NamedTemporaryFile() as tmp:
        ont = load()
        mapping = so.name_index(ont, tmp.name)
        assert mapping[so_id] == so_id
        assert mapping[name] == so_id


# 'transcript',
# 'pre-miRNA',
# 'hammerhead_ribozyme',
# 'riboswitch',
# 'ncRNA_gene',
# 'ribozyme',
# 'regulatory_region',
# 'primary_transcript',
# 'group_II_intron',
# 'five_prime_UTR',
# 'group_I_intron',
# 'cRISPR',
# 'attenuator',
# 'internal_ribosome_entry_site',
# 'autocatalytically_spliced_intron',
# 'three_prime_UTR',
# 'sECIS_element',
# 'repeat_unit',
# 'antisense',
# 'mRNA_region',
# 'mature_transcript',
# 'cis_regulatory_frameshift_element',
# 'recoding_stimulatory_region',
# 'rNA_sequence_secondary_structure',
# 'rRNA_gene',
