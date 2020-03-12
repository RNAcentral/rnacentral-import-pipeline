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


import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.psi_mi import tab


@pytest.mark.parametrize('raw,count', [
    ('psi-mi:q94533_drome(display_long)', 1),
    ('psi-mi:q94533_drome(display_long)|uniprotkb:S6k(gene name)', 2),
    ('psi-mi:q94533_drome(display_long)|uniprotkb:S6k(gene name)|psi-mi:S6k(display_short)|uniprotkb:S6k-RA(gene name synonym)|uniprotkb:10539(gene name synonym)|uniprotkb:7-10(gene name synonym)|uniprotkb:7084(gene name synonym)|uniprotkb:DS6K(gene name synonym)|uniprotkb:Dmel\CG10539(gene name synonym)|uniprotkb:Dp70S6k(gene name synonym)|uniprotkb:Dp70[s6k](gene name synonym)|uniprotkb:Dp70s6k(gene name synonym)|uniprotkb:Pk64F(gene name synonym)|uniprotkb:Pk?6(gene name synonym)|uniprotkb:S6K(gene name synonym)|uniprotkb:S6K1(gene name synonym)|uniprotkb:"S6k|S6K"(gene name synonym)|uniprotkb:dS6K(gene name synonym)|uniprotkb:dS6k(gene name synonym)|uniprotkb:dp70[S6k](gene name synonym)|uniprotkb:dp70s6k(gene name synonym)|uniprotkb:dps6k(gene name synonym)|uniprotkb:ds6k(gene name synonym)|uniprotkb:"fs(3)07084"(gene name synonym)|uniprotkb:"l(3)07084"(gene name synonym)|uniprotkb:p-S6K(gene name synonym)|uniprotkb:p70(gene name synonym)|uniprotkb:p70 S6K(gene name synonym)|uniprotkb:p70S6K(gene name synonym)|uniprotkb:p70[S6 kinase](gene name synonym)|uniprotkb:p70[S6K](gene name synonym)|uniprotkb:p70[S6k](gene name synonym)|uniprotkb:p70[S6kinase](gene name synonym)|uniprotkb:p70s6K(gene name synonym)|uniprotkb:s6k(gene name synonym)|uniprotkb:s6k11(gene name synonym)|uniprotkb:Dmel_CG10539(orf name)|uniprotkb:CG10539(orf name)|uniprotkb:d-S6K(gene name synonym)|uniprotkb:dP70S6K(gene name synonym)|uniprotkb:dS6K1(gene name synonym)', 41),
    ('"mimix:comment":not sure if the protein origin is right determined, my suggestion is based on fig.1c and 1d description|comment:"DIP protein Q05769 original sequence version: 98"', 2),
])
def test_parses_complex_field(raw, count):
    value = tab.identifiers(raw)
    assert len(value) == count


@pytest.mark.parametrize('raw,expected', [
    ('psi-mi:"MI:0396"(predetermined participant)', [
        data.InteractionIdentifier('psi-mi', 'MI:0396', 'predetermined participant')
    ]),
    ('uniprotkb:O43426', [
        data.InteractionIdentifier('uniprotkb', 'O43426', None)
    ]),
    ('uniprotkb:P34708-1', [
        data.InteractionIdentifier('uniprotkb', 'P34708-1', None)
    ]),
    ('go:"GO:0005525"(GTP binding)', [
        data.InteractionIdentifier('go', "GO:0005525", 'GTP binding')
    ]),
    ('rcsb pdb:1DYN', [
        data.InteractionIdentifier('rcsb pdb', '1DYN', None)
    ]),
    ('interpro:IPR003017(Amphiphysin, isoform 1)', [
        data.InteractionIdentifier('interpro', 'IPR003017', 'Amphiphysin, isoform 1')
    ]),
    ('refseq:NP_647477.1', [
        data.InteractionIdentifier('refseq', 'NP_647477.1', None)
    ]),
    ('go:"GO:0000977"(RNA polymerase II regulatory region sequence-specific DNA binding)', [
        data.InteractionIdentifier('go', 'GO:0000977', 'RNA polymerase II regulatory region sequence-specific DNA binding')
    ]),
    ('uniprotkb:"S6k|S6K"(gene name synonym)', [
        data.InteractionIdentifier('uniprotkb', 'S6k|S6K', 'gene name synonym')
    ]),
    ('psi-mi:"MI:0000"(a cv term)', [
        data.InteractionIdentifier('psi-mi', 'MI:0000', 'a cv term')
    ]),
    ('psi-mi:"MI:0000"("I can now use braces ()()() or pipes ||| here and ::colons::")', [
        data.InteractionIdentifier('psi-mi', 'MI:0000', 'I can now use braces ()()() or pipes ||| here and ::colons::')
    ]),
    (r'uniprotkb:P12345("a \"nice\" protein")', [
        data.InteractionIdentifier('uniprotkb', 'P12345', 'a "nice" protein')
    ]),
])
def test_parses_data_correctly(raw, expected):
    assert tab.identifiers(raw) == expected


@pytest.mark.parametrize('filename,count', [
    ('data/intact/sample.txt', 106),
    ('data/intact/problems.txt', 6),
    ('data/intact/quoting-issue.txt', 2),
])
def test_can_parse_all_data(filename, count):
    with open(filename, 'r') as raw:
        assert sum(1 for x in tab.parse(raw)) == count
