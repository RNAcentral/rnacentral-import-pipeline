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

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tasks.export.ftp.fasta.nhmmer import NHmmerIncludedExport
from tasks.export.ftp.fasta.nhmmer import NHmmerExcludedExport


@pytest.mark.parametrize('sequence,accepted', [
    ('', False),
    ('A', True),
    ('A-', False),
    ('ABCDGHKMNRSTVWXYU', True),
    ('GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAF', False),
    ('GCCCGGAUAGCUCAGUCGGUAGAGCAGGGGAUUGAAAAUCCCCGUGUCCUUGGUUCGAUUCCGAGUCCGGGCACCAF', False),
    ('GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCAFC', False),
    ('AUUIGAAAUC', False),
    ('CGGUGAIAAGGG', False),
    ('UUCCUCGUGGCCCAAUGGUCACGGCGUCUGGCUICGAACCAGAAGAUUCCAGGUUCAAGUCCUGGCGGGGAAGCCA', False),
    ('GGCUICGAACC', False),
    ('AUUIGAAAUCU', False),
    ('CUCGGCUICGAACCGAG', False),
    ('CUCGGXUICGAACCGAG', False),
])
def test_hmmer_class_make_correct_decisions(sequence, accepted):
    record = SeqRecord(Seq(sequence))
    included = NHmmerIncludedExport()
    excluded = NHmmerExcludedExport()
    assert included.is_valid_record(record) is accepted
    assert excluded.is_valid_record(record) is (not accepted)
