# -*- coding: utf-8 -*-

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CODING = (ROOT / "workflows" / "databases" / "ensembl_mrna.nf").read_text()
ENSEMBL = (ROOT / "workflows" / "databases" / "ensembl.nf").read_text()


def process(text, name):
    return text.split(f"process {name} {{", 1)[1].split("\n}\n", 1)[0]


def tuple_arity(text, name, section):
    body = process(text, name).split(f"{section}:", 1)[1]
    line = next(l for l in body.splitlines() if "tuple" in l)
    return len(re.findall(r"\b(?:val|path)\(", line))


# Unchanged genesets are skipped only if the diff reads the coding database's own signatures.
def test_coding_genesets_are_diffed_against_their_own_signatures():
    assert "--database ENSEMBL_CODING" in process(CODING, "coding_source_diff")
    assert "--database ENSEMBL_CODING" in process(CODING, "coding_delta_diff")
    assert "--source" in process(CODING, "parse_mrna")


# The loader only goes DELTA, and retires dropped accessions, if these files reach it.
def test_the_delta_files_reach_the_loader():
    tail = CODING.split("workflow ensembl_mrna {", 1)[1]
    for out in ["coding_delta_diff.out.manifest", "coding_delta_diff.out.deletions"]:
        assert out in tail
    assert "coding_source_diff.out.sources" in tail


# fetch_species_data grew to six fields on dev and parse_mrna still took three.
def test_parse_mrna_takes_what_fetch_species_data_emits():
    emitted = tuple_arity(ENSEMBL, "fetch_species_data", "output")
    assert tuple_arity(CODING, "parse_mrna", "input") == emitted


# Column 4 is the assembly id; a species with no homology file would wget an empty url.
def test_homology_is_fetched_only_where_there_is_a_url():
    assert "filter { row -> row[5] }" in CODING


# The EBI FTP refuses parallel downloads; the homology files are mirrored on /nfs too.
def test_homology_is_copied_from_nfs_on_a_datamover():
    body = process(CODING, "fetch_homology")
    assert "queue 'datamover'" in body
    assert "params.databases.ensembl.nfs_path" in body
    assert "wget" not in body
