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

from databases.pdb import helpers


@pytest.mark.parametrize('product,expected', [
    ('SRP RNA DOMAIN IV', 'SRP_RNA'),
    ('PRE-TRNA BULGE-HELIX-BULGE MOTIF', 'tRNA'),
    ('tmRNA (63-MER)', 'tmRNA'),
    ('U2 snRNA', 'snRNA'),
    ('HAMMERHEAD RIBOZYME-RNA STRAND', 'hammerhead_ribozyme'),
    ('RNA (HEPATITIS DELTA VIRUS GENOMIC RIBOZYME)', 'ribozyme'),
    ('U65 H/ACA snoRNA', 'ncRNA'),
    ("RNA (5'-R(*CP*GP*CP*GP*AP*AP*UP*UP*AP*GP*CP*G)-3')", 'misc_RNA'),
    ('5S RRNA LOOP D/LOOP E', 'rRNA'),
    ('5S RIBOSOMAL RNA', 'rRNA'),
    ('5.8S RRNA', 'rRNA'),
    ("16S RRNA (5'-R(*GP*GP*CP*GP*UP*CP*AP*CP*AP*CP*CP*UP*UP*C)-3')", 'rRNA'),
    ('RNA (16S RNA)', 'rRNA'),
    ('Helix 39 of 18S rRNA', 'rRNA'),
    ('18S ribosomal RNA', 'rRNA'),
    ('23S RRNA', 'rRNA'),
    ('HELIX 95 OF 23S RRNA', 'rRNA'),
    ('23S RRNA FRAGMENT', 'rRNA'),
    ('28S rRNA', 'rRNA'),
    ('sarcin/ricin 28S rRNA', 'rRNA'),
    ('30S 16S ribosomal RNA', 'rRNA'),
    ('30S RNA helix 8', 'rRNA'),
    ('40S ribosomal RNA fragment', 'rRNA'),
    ('40S ribosomal RNA', 'rRNA'),
    ('40S WHEAT GERM RIBOSOME protein 4', 'rRNA'),
    ('60S ribosomal RNA fragment', 'rRNA'),
    ('60S ribosomal RNA', 'rRNA'),
])
def test_can_compute_correct_rna_types(product, expected):
    assert helpers.rna_type({'compound': product}) == expected


def test_gets_given_taxid():
    assert(helpers.taxid({'taxonomyId': '562'})) == 562


def test_uses_synthenic_if_given_no_taxid():
    assert(helpers.taxid({'taxonomyId': ''})) == 32630


def test_can_detect_if_sequence_is_not_rna():
    data = {
        'chainLength': '25',
        'chainId': 'U',
        'db_name': '',
        'classification': 'RIBOSOME',
        'sequence': 'GKGDRRTRRGKIWRGTYGKYRPRKK',
        'resolution': '3.3',
        'structureId': '5WNP',
        'structureTitle': 'Crystal Structure of 30S ribosomal subunit from Thermus thermophilus',
        'emdbId': '',
        'releaseDate': '2018-02-21',
        'source': 'Thermus thermophilus',
        'ndbId': '',
        'db_id': '',
        'experimentalTechnique': 'X-RAY DIFFRACTION',
        'compound': "RNA (5'-R(*AP*AP*AP*UP*UP*U)-3')",
        'entityId': '21',
        'entityMacromoleculeType': 'Polyribonucleotide (RNA)',
        'taxonomyId': '274'
    }
    with pytest.raises(helpers.InvalidSequence):
        helpers.sequence(data)
