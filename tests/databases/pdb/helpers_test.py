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

from rnacentral_pipeline.databases.pdb.data import ChainInfo
from rnacentral_pipeline.databases.pdb import helpers
from rnacentral_pipeline.databases.pdb import fetch


def load(pdb_id: str, chain_id: str) -> ChainInfo:
    chains = fetch.rna_chains(pdb_ids=[pdb_id.lower()])
    return next(c for c in chains if c.chain_id == chain_id)


@pytest.mark.parametrize(
    "product,expected",
    [
        ("SRP RNA DOMAIN IV", "SRP_RNA"),
        ("PRE-TRNA BULGE-HELIX-BULGE MOTIF", "tRNA"),
        ("tmRNA (63-MER)", "tmRNA"),
        ("U2 snRNA", "snRNA"),
        ("HAMMERHEAD RIBOZYME-RNA STRAND", "hammerhead_ribozyme"),
        ("RNA (HEPATITIS DELTA VIRUS GENOMIC RIBOZYME)", "ribozyme"),
        ("U65 H/ACA snoRNA", "ncRNA"),
        ("RNA (5'-R(*CP*GP*CP*GP*AP*AP*UP*UP*AP*GP*CP*G)-3')", "misc_RNA"),
        ("5S RRNA LOOP D/LOOP E", "rRNA"),
        ("5S RIBOSOMAL RNA", "rRNA"),
        ("5.8S RRNA", "rRNA"),
        ("16S RRNA (5'-R(*GP*GP*CP*GP*UP*CP*AP*CP*AP*CP*CP*UP*UP*C)-3')", "rRNA"),
        ("RNA (16S RNA)", "rRNA"),
        ("Helix 39 of 18S rRNA", "rRNA"),
        ("18S ribosomal RNA", "rRNA"),
        ("23S RRNA", "rRNA"),
        ("HELIX 95 OF 23S RRNA", "rRNA"),
        ("23S RRNA FRAGMENT", "rRNA"),
        ("28S rRNA", "rRNA"),
        ("sarcin/ricin 28S rRNA", "rRNA"),
        ("30S 16S ribosomal RNA", "rRNA"),
        ("30S RNA helix 8", "rRNA"),
        ("40S ribosomal RNA fragment", "rRNA"),
        ("40S ribosomal RNA", "rRNA"),
        ("40S WHEAT GERM RIBOSOME protein 4", "rRNA"),
        ("60S ribosomal RNA fragment", "rRNA"),
        ("60S ribosomal RNA", "rRNA"),
    ],
)
def test_can_compute_correct_rna_types(product: str, expected):
    assert helpers.compound_rna_type(product) == expected


@pytest.mark.parametrize(
    "pdb,chain,expected",
    [
        ("7mky", "A", True),
        ("7lyj", "A", True),
        ("5U3G", "B", True),
        ("2L1V", "A", True),
        ("6VAR", "A", True),
        ("4Y1I", "A", True),
        ("4Y1I", "B", True),
        ("4Y1J", "A", True),
        ("4Y1J", "B", True),
        ("4Y1M", "A", True),
        ("4Y1M", "B", True),
        ("7MKY", "A", True),
        ("7LYJ", "A", True),
        ("7MLW", "F", True),
    ],
)
def test_can_detect_if_is_ncrna(pdb, chain, expected):
    info = load(pdb, chain)
    assert helpers.is_ncrna(info) == expected
