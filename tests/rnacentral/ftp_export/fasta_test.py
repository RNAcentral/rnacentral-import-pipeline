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

from rnacentral.export.ftp import fasta
from rnacentral.psql import PsqlWrapper

from tasks.config import db

from tests.tasks.helpers import count


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
def test_can_classify_sequences_for_nhmmer(sequence, accepted):
    record = SeqRecord(Seq(sequence))
    assert fasta.is_valid_nhmmer_record(record) is accepted


def test_can_produce_records_from_a_query():
    sql = "select upi as id, md5 as description, 'AAAA' as sequence from rna order by id asc limit 3"
    psql = PsqlWrapper(db())
    lines = [fasta.as_record(l) for l in psql.copy_to_iterable(sql)]
    # Biopython doesn't implement SeqRecord.__eq__, ugh
    assert len(lines) == 3
    assert lines[0].id == 'URS0000000001'
    assert lines[0].description == '6bba097c8c39ed9a0fdf02273ee1c79a'
    assert lines[0].seq == Seq('AAAA')

    assert lines[1].id == 'URS0000000002'
    assert lines[1].description == '1fe2f0e3c3a2d6d708ac98e9bfb1d7a8'
    assert lines[1].seq == Seq('AAAA')

    assert lines[2].id == 'URS0000000003'
    assert lines[2].description == '7bb11d0572ff6bb42427ce74450ba564'
    assert lines[2].seq == Seq('AAAA')


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_active_sql_produces_correct_counts():
    assert count(fasta.ACTIVE_SQL) == 11368718


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_active_produces_correct_number_of_records():
    assert sum(1 for _ in fasta.active(db())) == 11368718



# @pytest.mark.slowtest
@pytest.mark.skip()
def test_active_has_correct_first_record():
    data = next(fasta.active(db()))
    assert data.id == 'URS0000000001'
    assert data.description == 'environmental samples uncultured bacterium rRNA'
    assert data.seq == Seq('AUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGUAGAGAGAAGCUUGCUUCUCUUGAGAGCGGCGGACGGGUGAGUAAUGCCUAGGAAUCUGCCUGGUAGUGGGGGAUAACGCUCGGAAACGGACGCUAAUACCGCAUACGUCCUACGGGAGAAAGCAGGGGACCUUCGGGCCUUGCGCUAUCAGAUGAGC')


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_species_sql_produces_correct_counts():
    assert count(fasta.ACITVE_SPECIES_SQL) == 13788772


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_species_produces_correct_number_of_records():
    assert sum(1 for _ in fasta.species(db())) == 13788772


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_species_has_correct_first_record():
    data = next(fasta.species(db()))
    assert data.id == 'URS000050A1B4_77133'
    assert data.description == 'uncultured bacterium partial 16S ribosomal RNA'
    assert data.seq == Seq('TGAGGAATATTGGTCAATGGGCGAGAGCCTGAAACCAGCCAAGTAGCGTGAAGGAAGACTGCCCTATGGGTTGTAAACTTCTTTTATAAGGGAATAAAGAGCGCCACGTGTGGTGTGTTGTATGTACCTTATGAATAAGCATCGGCTAATTCCGTGCC')


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_inactive_sql_produces_correct_counts():
    assert count(fasta.INACTIVE_SQL) == 2079637


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_inactive_produces_correct_number_of_records():
    assert sum(1 for _ in fasta.inactive(db())) == 2079637


# @pytest.mark.slowtest
@pytest.mark.skip()
def test_inactive_has_correct_first_record():
    data = next(fasta.inactive(db()))
    assert data.id == 'URS00006E9C09'
    assert data.description == 'Saccoglossus kowalevskii rRNA'
    assert data.seq == Seq('GCCTACTGCCAAACCATTGCAAAATACACCAGTTCTCATCCGATCACTGAAGTTAAGCTGCATCGGGCGTGGTTAGTACTTGCATGGGAGACGGGTAGGGAATACCACGTGCAGTAGGCT')
