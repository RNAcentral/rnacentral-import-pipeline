# -*- coding: utf-8 -*-

"""
piRBase gives us a bare piRNA id and a sequence, so taxid, rna_type and
description are all inferred here and this pins the inference.

The filtering is the part worth guarding. The full sets are far too large to
import whole and are cut down to sequences RNAcentral already holds, while the
gold standard sets must bypass that filter entirely -- their ids barely overlap
the full sets, so a new gold sequence has nothing to match against. Getting
that backwards would silently drop almost all the curated data.
"""

from pathlib import Path

import pytest

from rnacentral_pipeline.databases.pirbase import parser
from rnacentral_pipeline.databases.pirbase.data import SPECIES, species, urls

CEL = Path("data/pirbase/cel.fa")
HSA_GOLD = Path("data/pirbase/hsa.gold.fa")

BASE = "http://bigdata.ibp.ac.cn/piRBase/download/v3.0"


def parse(path, code, known=None):
    with path.open("r") as handle:
        return list(parser.parse(handle, species(code), known_path=known))


def test_builds_entries_from_a_full_set():
    entries = parse(CEL, "cel")
    assert len(entries) == 10

    entry = entries[0]
    assert entry.accession == "piR-cel-1"
    assert entry.primary_id == "piR-cel-1"
    assert entry.ncbi_tax_id == 6239
    assert entry.database == "PIRBASE"
    assert entry.rna_type == "SO:0001035"
    assert entry.sequence == "TGGTACGTACGTTAACCGTGC"
    assert entry.description == "Caenorhabditis elegans piRNA piR-cel-1"
    assert entry.species == "Caenorhabditis elegans"
    # piRBase publishes BED only against assemblies RNAcentral no longer
    # carries, and maps these sequences itself.
    assert entry.regions == []
    assert entry.is_valid()


def test_gold_sets_are_not_filtered():
    entries = parse(HSA_GOLD, "hsa")
    assert len(entries) == 6
    assert entries[0].accession == "piR-hsa-9"
    assert entries[0].ncbi_tax_id == 9606


def test_full_sets_keep_only_known_md5s(tmp_path):
    md5s = tmp_path / "md5s"
    # piR-cel-2 and piR-cel-9, plus an md5 belonging to nothing in the file.
    md5s.write_text(
        "\n".join(
            [
                "d406fd080142ac1cff336f8aacabb20e",
                "35df1198a039dc4f22e2991101f520bb",
                "0" * 32,
            ]
        )
        + "\n"
    )
    known = tmp_path / "known.sqlite"
    parser.build_known(md5s, known)

    entries = parse(CEL, "cel", known=known)
    assert [e.accession for e in entries] == ["piR-cel-2", "piR-cel-9"]


def test_an_empty_known_index_keeps_nothing(tmp_path):
    md5s = tmp_path / "md5s"
    md5s.write_text("")
    known = tmp_path / "known.sqlite"
    parser.build_known(md5s, known)

    assert parse(CEL, "cel", known=known) == []


def test_covers_every_species_piRBase_publishes():
    assert len(SPECIES) == 43
    assert len({s.taxid for s in SPECIES.values()}) == 43


@pytest.mark.parametrize(
    "code,name,taxid",
    [
        # piRBase still uses the collection numbers these three had before they
        # were described; they were resolved through their assemblies, not by
        # name, and must not be "corrected" back to a Caenorhabditis sp.
        ("c26", "Caenorhabditis zanzibari", 2306312),
        ("c31", "Caenorhabditis uteleia", 2305860),
        ("c32", "Caenorhabditis sulstoni", 2305862),
        # A Plectus, despite the abbreviation looking like Panagrolaimus.
        ("psa", "Plectus sambesii", 2011161),
        # Macaca fascicularis, not mulatta; piRBase ships both.
        ("mfa", "Macaca fascicularis", 9541),
        ("mml", "Macaca mulatta", 9544),
    ],
)
def test_pins_the_taxids_that_needed_looking_up(code, name, taxid):
    assert species(code).name == name
    assert species(code).taxid == taxid


def test_generates_a_url_for_every_file():
    generated = list(urls(BASE))
    assert len(generated) == 49
    assert sum(1 for _, kind, _ in generated if kind == "gold") == 6

    by_kind = {(code, kind): url for code, kind, url in generated}
    # piRBase is inconsistent about whether the release is in the filename.
    assert by_kind[("hsa", "full")] == f"{BASE}/fasta/hsa.v3.0.fa.gz"
    assert by_kind[("rno", "full")] == f"{BASE}/fasta/rno.fa.gz"
    assert by_kind[("hsa", "gold")] == f"{BASE}/fasta/hsa.gold.fa.gz"


def test_rejects_an_unknown_species():
    with pytest.raises(ValueError):
        species("nope")
