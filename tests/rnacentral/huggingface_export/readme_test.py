# -*- coding: utf-8 -*-

from contextlib import contextmanager

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import huggingface


@pytest.fixture
def readme(tmp_path, monkeypatch):
    @contextmanager
    def fake_cursor(_url):
        yield None

    monkeypatch.setattr(huggingface, "cursor", fake_cursor)
    monkeypatch.setattr(huggingface, "get_active_sequence_count", lambda _cur: 42)
    specs = tmp_path / "specs.txt"
    specs.write_text("RNAcentral: 1")

    def build(*flags):
        result = CliRunner().invoke(
            huggingface.cli,
            ["create-readme", "27", str(specs), *flags],
            env={"PGDATABASE": "unused"},
        )
        assert result.exit_code == 0, result.output
        return (tmp_path / "README.md").read_text()

    monkeypatch.chdir(tmp_path)
    return build


def test_default_readme_declares_no_splits(readme):
    text = readme()
    assert "configs:" not in text
    assert "deduplicated" not in text


def test_deduplicated_split_is_declared_beside_train(readme):
    text = readme("--with-deduplicated")
    front_matter = text.split("---")[1]
    assert "split: train" in front_matter
    assert "path: active_sequences.parquet" in front_matter
    assert "split: deduplicated" in front_matter
    assert "path: deduplicated.parquet" in front_matter
    assert "canonical" in text
