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

from databases.ena import dr


def test_can_parse_a_dr_line():
    ans = dr.DBRef('miRBase', 'MI0031512', 'hsa-mir-3648-1')
    assert dr.parse_line('DR   miRBase; MI0031512; hsa-mir-3648-1.') == ans


def test_can_parse_a_dr_line_with_only_primary():
    ans = dr.DBRef('tmRNA-Website', 'Acary_marin_MBIC11', None)
    assert dr.parse_line('DR   tmRNA-Website; Acary_marin_MBIC11.') == ans


def test_can_extract_dr_lines():
    with open('data/ena/tpa/mirbase/entry.embl', 'rb') as raw:
        data = dr.mapping(raw)

    assert data == {
        'LM611181.1:1..180:precursor_RNA': [
            dr.DBRef('miRBase', 'MI0016048', 'hsa-mir-3648-1'),
            dr.DBRef('MD5', '08036e5a2a91e75299436501f4182050', None),
        ],
        'LM611890.1:1..180:precursor_RNA': [
            dr.DBRef('miRBase', 'MI0031512', 'hsa-mir-3648-2'),
            dr.DBRef('MD5', '08036e5a2a91e75299436501f4182050', None),
        ]
    }
