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

from rnacentral_pipeline.databases.rfam import utils


@pytest.mark.parametrize(
    "name,expected",
    [
        ("Cis52_sRNA", "other"),
        ("CbSR2", "antisense_RNA"),
        ("ppoRNA", "rRNA"),
        ("SNOR75", "snoRNA"),
        ("URE2_IRES", "other"),
        ("mir-393", "precursor_RNA"),
        ("GlsR20", "snoRNA"),
        ("Hammerhead_HH10", "hammerhead_ribozyme"),
        ("TTC28-AS1_3", "lncRNA"),
        ("ROSE_2", "other"),
        ("astro_FSE", "other"),
        ("RUF20", "other"),
        ("lactis-plasmid", "other"),
        ("leu-phe_leader", "other"),
        ("CRISPR-DR2", "other"),
        ("Metazoa_SRP", "SRP_RNA"),
        ("TP53TG1_1", "lncRNA"),
        ("group-II-D1D4-4", "autocatalytically_spliced_intron"),
        ("Vault", "vault_RNA"),
        ("tmRNA", "tmRNA"),
        ("CDKN2B-AS", "lncRNA"),
        ("7SK", "snRNA"),
        ("isrP", "other"),
        ("U2", "snRNA"),
    ],
)
def test_can_fetch_a_mapping_from_name_to_isndc(name, expected):
    with open("data/rfam/families.tsv") as raw:
        mapping = utils.name_to_insdc_type(raw)
    assert mapping[name] == expected


@pytest.mark.parametrize(
    "family_id,expected",
    [
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
        ("RF01414", "other"),
        ("RF01502", "SRP_RNA"),
        ("RF01577", "RNase_P_RNA"),
        ("RF01695", "antisense_RNA"),
        ("RF01807", "autocatalytically_spliced_intron"),
        ("RF01854", "SRP_RNA"),
        ("RF01855", "SRP_RNA"),
        ("RF01856", "SRP_RNA"),
        ("RF01857", "SRP_RNA"),
        ("RF01960", "rRNA"),
        ("RF02144", "other"),
        ("RF02276", "hammerhead_ribozyme"),
        ("RF02277", "hammerhead_ribozyme"),
        ("RF02357", "RNase_P_RNA"),
        ("RF02554", "rRNA"),
        ("RF02558", "antisense_RNA"),
        ("RF02647", "other"),
    ],
)
def test_can_fetch_a_mapping_from_id_to_isndc(family_id, expected):
    with open("data/rfam/families.tsv") as raw:
        mapping = utils.id_to_insdc_type(raw)
    assert mapping[family_id] == expected


def test_it_maps_a_known_lncRNA():
    with open("data/rfam/families.tsv") as raw:
        mapping = utils.id_to_insdc_type(raw)
    assert mapping["RF01800"] == "lncRNA"
