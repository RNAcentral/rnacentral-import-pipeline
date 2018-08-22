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

import pytest


@pytest.mark.paramterize('rna_id,expected', [
    ('URS00000C777C_1231464', "Alkalibacterium sp. 3.5P*23 small subunit ribosomal RNA (16S)"),
    ('URS000080DD59_32630', "5'-R(*UP*GP*(CBV)P*(CBV)P*AP*GP*UP*UP*CP*GP*CP*UP*GP*GP*C)-3' from (PDB 1QBP, chain E)"),
    ('URS000080DE07_32630', "RNA (5'-R(P*CP*GP*AP*UP*CP*GP*GP*GP*UP*GP*UP*C)-3') from (PDB 1RY1, chain Q)"),
    ('URS000080DE5B_32630', "RNA (5'-R(P*AP*UP*CP*GP*CP*GP*CP*CP*UP*GP*UP*G)-3') from (PDB 1RY1, chain R)"),
    ('URS000080E209_32630', "RNA (5'-R(*CP*CP*GP*CP*CP*GP*CP*GP*CP*CP*AP*(5BU)P*GP*CP*CP*UP*GP*UP*GP*GP*CP*GP... from (PDB 3MEI, chain B)"),
    ('URS000080E22B_308052', "tRNA (5'-D(*AP*UP*CP*CP*CP*CP*GP*UP*GP*UP*CP*CP*UP*UP*GP*GP*UP*UP*CP*G)-3') from Mitsuaria sp. 67 (PDB 4WT8, chain D2)"),
    ('URS000080E230_32630', "5'-D(*CP*AP*GP*CP*TP*AP*CP*TP*TP*GP*AP*GP*CP*T)-3' from (PDB 3H3V, chain P)"),
    ('URS0000A77852_32630', "RNA (5'-R(*(LCC)P*(LCC)P*(LCA)P*(LCG)P*AP*CP*UP*UP*AP*AP*GP*UP*CP*U)-3') from synthetic construct (PDB 5V0J, chain B)"),
    ('URS0000BC4693_9606', "RNA (5'-R(*GP*(CBV)P*CP*GP*GP*(6MZ)P*UP*GP*GP*C)-3') from Homo sapiens (PDB 5LR4, chain D)"),
    ('URS000080E10B_32630', "5'-R(*GP*CP*CP*GP*AP*AP*GP*CP*CP*(P5P)-3' from (PDB 1XV0, chain B)"),
    ('URS000080E135_274', "messenger RNA (5'-R(*AP*AP*UP*GP*UP*AP*G)-3') from Thermus thermophilus (PDB 4V9N, chain CV)"),
    ('URS000041AF00_274', 'RNA (77-MER) from Thermus thermophilus (PDB 4V7J, chain Bw)'),
])
def test_computes_good_masked_description(rna_id, expected):
    assert False
