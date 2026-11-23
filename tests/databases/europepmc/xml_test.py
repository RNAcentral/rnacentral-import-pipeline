# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import attr

from rnacentral_pipeline.databases.data import Reference
from rnacentral_pipeline.databases.europepmc import xml


def test_can_parse_xml_data_correctly():
    with open("data/publications/example.xml", "rb") as raw:
        data = list(xml.parse(raw))
        assert len(data) == 1
        assert attr.asdict(data[0]) == attr.asdict(
            Reference(
                authors="Xu Z, Han Y, Liu J, Jiang F, Hu H, Wang Y, Liu Q, Gong Y, Li X.",
                location="Scientific reports 5:12276 (2015)",
                title=(
                    "MiR-135b-5p and MiR-499a-3p Promote Cell "
                    "Proliferation and Migration in Atherosclerosis by Directly "
                    "Targeting MEF2C"
                ),
                pmid=26184978,
                doi="10.1038/srep12276",
                pmcid="PMC4505325",
            )
        )
