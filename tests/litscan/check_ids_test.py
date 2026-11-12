# -*- coding: utf-8 -*-

import importlib.util
import sys
import types
from pathlib import Path

from click.testing import CliRunner

ROOT = Path(__file__).resolve().parents[2]


def load_checker(monkeypatch):
    # The nltk word list is only installed in the container.
    words = types.SimpleNamespace(words=lambda: [])
    nltk = types.SimpleNamespace(corpus=types.SimpleNamespace(words=words))
    monkeypatch.setitem(sys.modules, "nltk", nltk)
    path = ROOT / "bin" / "litscan-check-ids.py"
    spec = importlib.util.spec_from_file_location("litscan_check_ids", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Papers cite Ensembl ids without their version, so a versioned search finds nothing.
def test_coding_ids_are_searched_without_their_version(monkeypatch, tmp_path):
    checker = load_checker(monkeypatch)
    raw = tmp_path / "results"
    raw.write_text("URS0000ABCD|9606|ENST00000269305.9|ENSG00000141510.18\n")
    out = tmp_path / "output"
    result = CliRunner().invoke(checker.main, ["ensembl_coding", str(raw), str(out)])
    assert result.exit_code == 0, result.output
    assert "ENST00000269305.9" not in out.read_text()
    assert "ENSG00000141510|ENST00000269305|URS0000ABCD_9606" in out.read_text()


def test_the_coding_query_searches_mrnas_by_their_ids_only():
    query = (
        ROOT / "workflows" / "litscan" / "queries" / "ensembl_coding.sql"
    ).read_text()
    assert "'ENSEMBL_CODING'" in query
    assert "feature_name = 'mRNA'" in query
    assert "locus_tag" not in query
