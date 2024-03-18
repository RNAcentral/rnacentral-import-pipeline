# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
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

from pathlib import Path, PosixPath

import pytest

from rnacentral_pipeline.databases.expressionatlas import configuration


@pytest.fixture(scope="module")
def data():
    with open("data/expressionatlas/E-CURD-113-configuration.xml", "r") as raw:
        data = configuration.parse(raw)
    print(data)
    return data


@pytest.fixture(scope="module")
def path():
    return "data/expressionatlas/E-CURD-113-configuration.xml.xml"


def test_can_parse_ea_config(data):

    assert data["experimentType"] == "rnaseq_mrna_differential"
    assert len(data["analytics"]["assay_groups"]) == 4


def test_can_extract_exp_name(path):
    exp_name = configuration.extract_experiment_name(path)
    assert exp_name == "E-CURD-113"


@pytest.mark.parametrize(
    "paths,expected",
    [
        (
            [
                PosixPath("data/expressionatlas/E-ATMX-6-configuration.xml"),
                PosixPath("data/expressionatlas/E-CURD-113-configuration.xml"),
                PosixPath("data/expressionatlas/E-MTAB-3103-configuration.xml"),
            ],
            {
                "E-ATMX-6": {
                    "experimentType": "microarray_1colour_mrna_differential",
                    "r_data": "1",
                    "analytics": {
                        "array_design": "A-AFFY-2",
                        "assay_groups": {
                            "g1": ["H_WT-3", "H_WT-2", "H_WT-1"],
                            "g2": ["H_Myb28KO-1", "H_Myb28KO-3", "H_Myb28KO-2"],
                            "g3": ["H_Myb29KO-3", "H_Myb29KO-1", "H_Myb29KO-2"],
                        },
                        "contrasts": {
                            "g1_g2": {
                                "name": "'MYB28 knock out' vs 'wild type'",
                                "reference_assay_group": "g1",
                                "test_assay_group": "g2",
                            },
                            "g1_g3": {
                                "name": "'MYB29 knock out' vs 'wild type'",
                                "reference_assay_group": "g1",
                                "test_assay_group": "g3",
                            },
                        },
                    },
                },
                "E-CURD-113": {
                    "experimentType": "rnaseq_mrna_differential",
                    "r_data": "1",
                    "analytics": {
                        "array_design": None,
                        "assay_groups": {
                            "g1": [
                                "SRR1041550",
                                "SRR1041551",
                                "SRR1041548",
                                "SRR1041549",
                            ],
                            "g2": [
                                "SRR1041552",
                                "SRR1041553",
                                "SRR1041555",
                                "SRR1041554",
                            ],
                            "g3": ["SRR1041556", "SRR1041558", "SRR1041557"],
                            "g4": ["SRR1041559", "SRR1041561", "SRR1041560"],
                        },
                        "contrasts": {
                            "g1_g2": {
                                "name": "'fzt mutation' vs 'control' in 'seedling development stage; whole organism'",
                                "reference_assay_group": "g1",
                                "test_assay_group": "g2",
                            },
                            "g3_g4": {
                                "name": "'fzt mutation' vs 'control' in 'sporophyte vegetative stage; tassel primordia'",
                                "reference_assay_group": "g3",
                                "test_assay_group": "g4",
                            },
                        },
                    },
                },
                "E-MTAB-3103": {
                    "experimentType": "rnaseq_mrna_baseline",
                    "r_data": "1",
                    "analytics": {
                        "array_design": None,
                        "assay_groups": {
                            "g1": ["ERR687868", "ERR687864", "ERR687867"],
                            "g2": ["ERR687863", "ERR687866", "ERR687861"],
                            "g3": ["ERR687862", "ERR687869", "ERR687865"],
                        },
                    },
                },
            },
        )
    ],
)
def test_can_build_lookup(paths, expected):
    lookup_dict = configuration.build_exp_config_lookup(paths)

    print(lookup_dict)
    assert lookup_dict.keys() == expected.keys()
    assert lookup_dict["E-CURD-113"] == expected["E-CURD-113"]
    assert lookup_dict["E-ATMX-6"] == expected["E-ATMX-6"]
    assert lookup_dict["E-MTAB-3103"] == expected["E-MTAB-3103"]
