# -*- coding: utf-8 -*-

"""
RiboVision was imported by scraping rnacentral_mapping.html, keyed on PDB
chains. apollo.chemistry.gatech.edu is gone and RiboVision2 serves a JavaScript
app with no equivalent page, so the models now come from the copy R2DT ships.

These cover the part that needs no taxonomy service: reading bpseq2fasta output.
Entries carry no organelle, since nothing stored one before.
"""

import pytest

from rnacentral_pipeline.databases.helpers import r2dt
from rnacentral_pipeline.databases.ribovision import helpers

LSU = """> generated from /rna/auto-traveler/data/ribovision/bpseq/EC_LSU_3D.bpseq by bpseq2fasta
GGUUAAGCGACUAAGCGUACACGGUGGAUGCC
..((((((......((((((((((.....(((
"""

SSU = """> generated from /rna/auto-traveler/data/ribovision/bpseq/HS_SSU_3D.bpseq by bpseq2fasta
UACCUGGUUGAUCCUGCCAGUAGCAUAUGCUU
.....((((((...((((((((((((....((
"""


@pytest.fixture()
def r2dt_data(tmp_path):
    for subunit, content in (("lsu", LSU), ("ssu", SSU)):
        directory = tmp_path / f"ribovision-{subunit}" / "bpseq"
        directory.mkdir(parents=True)
        name = "EC_LSU_3D" if subunit == "lsu" else "HS_SSU_3D"
        (directory / f"{name}.fasta").write_text(content)
    return tmp_path


def test_reads_both_subunits(r2dt_data):
    found = {record.id: str(record.seq) for record in helpers.fasta_entries(r2dt_data)}
    assert found == {
        "EC_LSU_3D": "GGUUAAGCGACUAAGCGUACACGGUGGAUGCC",
        "HS_SSU_3D": "UACCUGGUUGAUCCUGCCAGUAGCAUAUGCUU",
    }


def test_drops_the_structure_line(r2dt_data):
    for record in helpers.fasta_entries(r2dt_data):
        assert "(" not in str(record.seq)


def test_rejects_a_header_without_a_model_id(tmp_path):
    directory = tmp_path / "bpseq"
    directory.mkdir(parents=True)
    (directory / "broken.fasta").write_text("> no model here\nACGU\n....\n")
    with pytest.raises(ValueError):
        list(r2dt.fasta_entries([directory]))


def test_tolerates_a_missing_trailing_newline(tmp_path):
    """The old reader unpacked exactly three lines and died without one."""
    directory = tmp_path / "bpseq"
    directory.mkdir(parents=True)
    (directory / "x.fasta").write_text(
        "> generated from a/b/EC_LSU_3D.bpseq by bpseq2fasta\nACGU\n...."
    )
    assert [r.id for r in r2dt.fasta_entries([directory])] == ["EC_LSU_3D"]
