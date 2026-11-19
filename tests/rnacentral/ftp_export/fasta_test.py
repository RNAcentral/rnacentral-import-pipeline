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

import io

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rnacentral_pipeline.rnacentral.ftp_export import fasta


@pytest.mark.parametrize(
    "sequence,accepted",
    [
        ("", False),
        ("A", True),
        ("A-", False),
        ("ABCDGHKMNRSTVWXYU", True),
        (
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAF",
            False,
        ),
        (
            "GCCCGGAUAGCUCAGUCGGUAGAGCAGGGGAUUGAAAAUCCCCGUGUCCUUGGUUCGAUUCCGAGUCCGGGCACCAF",
            False,
        ),
        (
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAFC",
            False,
        ),
        ("AUUIGAAAUC", False),
        ("CGGUGAIAAGGG", False),
        (
            "UUCCUCGUGGCCCAAUGGUCACGGCGUCUGGCUICGAACCAGAAGAUUCCAGGUUCAAGUCCUGGCGGGGAAGCCA",
            False,
        ),
        ("GGCUICGAACC", False),
        ("AUUIGAAAUCU", False),
        ("CUCGGCUICGAACCGAG", False),
        ("CUCGGXUICGAACCGAG", False),
    ],
)
def test_can_classify_sequences_for_nhmmer(sequence, accepted):
    record = SeqRecord(Seq(sequence))
    assert fasta.is_valid_nhmmer_record(record) is accepted


FASTA = ">good\nACGU\n>bad\nACGI\n>also_good\nNNNN\n"


def test_valid_nhmmer_keeps_only_the_accepted_records():
    out = io.StringIO()
    fasta.valid_nhmmer(io.StringIO(FASTA), out)
    assert [r.id for r in SeqIO.parse(io.StringIO(out.getvalue()), "fasta")] == [
        "good",
        "also_good",
    ]


def test_invalid_nhmmer_keeps_only_the_rejected_records():
    out = io.StringIO()
    fasta.invalid_nhmmer(io.StringIO(FASTA), out)
    assert [r.id for r in SeqIO.parse(io.StringIO(out.getvalue()), "fasta")] == ["bad"]
