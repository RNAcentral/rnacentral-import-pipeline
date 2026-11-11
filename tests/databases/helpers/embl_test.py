# -*- coding: utf-8 -*-

from rnacentral_pipeline.databases.helpers import embl

EMBL = """\
FH   Key             Location/Qualifiers
FT   gene            1..30
FT                   /gene="G1"
FT   CDS             join(1..10,
FT                   20..30)
FT                   /translation="MSGESA"
FT   misc_RNA        1..30
FT                   /note="snoRNA"
FT   CDS             complement(1..9)
XX
SQ   Sequence 30 BP;
     acgt
//
"""


def stripped(text):
    return "".join(embl.without_cds(text.splitlines(keepends=True)))


def test_it_drops_a_cds_with_all_its_continuation_lines():
    out = stripped(EMBL)
    assert "CDS" not in out
    assert "MSGESA" not in out
    assert "20..30" not in out


def test_it_keeps_the_features_around_a_cds():
    out = stripped(EMBL)
    assert 'FT   gene            1..30\nFT                   /gene="G1"\n' in out
    assert 'FT   misc_RNA        1..30\nFT                   /note="snoRNA"\n' in out


def test_it_keeps_the_sequence_after_a_final_cds():
    assert stripped(EMBL).endswith("XX\nSQ   Sequence 30 BP;\n     acgt\n//\n")


def test_it_leaves_a_file_without_cds_alone():
    text = EMBL.replace("CDS            ", "mRNA           ")
    assert stripped(text) == text
