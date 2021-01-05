# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.data import RnaType
from rnacentral_pipeline.databases.sequence_ontology import tree


@pytest.fixture(scope="module")
def so_tree():
    return tree.load_ontology(tree.REMOTE_ONTOLOGY)


@pytest.mark.parametrize(
    "name,so_id",
    [
        ("RNase_MRP_RNA", "SO:0000385"),
        ("RNase_P_RNA", "SO:0000386"),
        ("SRP_RNA", "SO:0000590"),
        ("Y_RNA", "SO:0000405"),
        ("antisense_RNA", "SO:0000644"),
        ("autocatalytically_spliced_intron", "SO:0000588"),
        ("guide_RNA", "SO:0000602"),
        ("hammerhead_ribozyme", "SO:0000380"),
        ("lncRNA", "SO:0001877"),
        ("miRNA", "SO:0000276"),
        ("ncRNA", "SO:0000655"),
        ("misc_RNA", "SO:0000673"),
        ("other", "SO:0000655"),
        ("precursor_RNA", "SO:0000185"),
        ("piRNA", "SO:0001035"),
        ("rasiRNA", "SO:0000454"),
        ("ribozyme", "SO:0000374"),
        ("scRNA", "SO:0000013"),
        ("scaRNA", "SO:0002095"),
        ("siRNA", "SO:0000646"),
        ("snRNA", "SO:0000274"),
        ("snoRNA", "SO:0000275"),
        ("telomerase_RNA", "SO:0000390"),
        ("tmRNA", "SO:0000584"),
        ("vault_RNA", "SO:0000404"),
        ("rRNA", "SO:0000252"),
        ("tRNA", "SO:0000253"),
        ("pre_miRNA", "SO:0001244"),
        ("sRNA", "SO:0000655"),
    ],
)
def test_can_build_from_insdc_name(so_tree, name, so_id):
    assert RnaType.from_insdc_term(so_tree, name).so_term.so_id == so_id


@pytest.mark.parametrize("so_id", [
    ("SO:0000252"),
    ("SO:0001244"),
    ("SO:0000650"),
    ("SO:0001000"),
])
def test_can_build_from_so_id(so_tree, so_id):
    assert RnaType.from_so_id(so_tree, so_id).so_term.so_id == so_id


@pytest.mark.parametrize("name,so_id", [
    ("rRNA", "SO:0000252"),
    ("pre_miRNA", "SO:0001244"),
    ("small_subunit_rRNA", "SO:0000650"),
    ("rRNA_16S", "SO:0001000"),
])
def test_can_build_from_so_name(so_tree, name, so_id):
    assert RnaType.from_so_id(so_tree, name).so_term.so_id == so_id
