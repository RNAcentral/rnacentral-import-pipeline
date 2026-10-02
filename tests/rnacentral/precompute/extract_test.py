# -*- coding: utf-8 -*-

from unittest.mock import patch

from click.testing import CliRunner

from rnacentral_pipeline.cli.precompute import cli as precompute_cli
from rnacentral_pipeline.rnacentral.precompute import extract


def test_substitute_range_replaces_min_max_with_null_when_unranged():
    sql = "SELECT * FROM t WHERE (:min IS NULL OR id BETWEEN :min AND :max)"

    result = extract._substitute_range(sql, None)

    assert result == "SELECT * FROM t WHERE (NULL IS NULL OR id BETWEEN NULL AND NULL)"


def test_substitute_range_replaces_min_max_with_real_integers():
    sql = "SELECT * FROM t WHERE (:min IS NULL OR id BETWEEN :min AND :max)"

    result = extract._substitute_range(sql, (10, 20))

    assert result == "SELECT * FROM t WHERE (10 IS NULL OR id BETWEEN 10 AND 20)"


def test_substitute_range_handles_multiple_occurrences_of_each_placeholder():
    # r2dt-hits.sql-shaped case: :min/:max each appear more than once
    # (once in the nullable check, once in the BETWEEN) - a naive
    # single-replace would miss the second occurrence.
    sql = ":min :min :max :max"

    result = extract._substitute_range(sql, (1, 2))

    assert result == "1 1 2 2"


def test_extract_query_cli_wires_through_to_extract_query(tmp_path, monkeypatch):
    monkeypatch.setenv("PGDATABASE", "postgresql://user:pass@host/db")
    sql_path = tmp_path / "basic.sql"
    sql_path.write_text("SELECT 1")
    output_path = tmp_path / "basic.parquet"

    with patch(
        "rnacentral_pipeline.cli.precompute.pre_extract.extract_query"
    ) as mock_extract:
        result = CliRunner().invoke(
            precompute_cli,
            ["extract-query", str(sql_path), str(output_path), "--range", "1", "100"],
        )

    assert result.exit_code == 0, result.output
    mock_extract.assert_called_once_with(
        sql_path, output_path, "postgresql://user:pass@host/db", range=(1, 100)
    )


def test_extract_query_cli_without_range_passes_none(tmp_path, monkeypatch):
    monkeypatch.setenv("PGDATABASE", "postgresql://user:pass@host/db")
    sql_path = tmp_path / "basic.sql"
    sql_path.write_text("SELECT 1")
    output_path = tmp_path / "basic.parquet"

    with patch(
        "rnacentral_pipeline.cli.precompute.pre_extract.extract_query"
    ) as mock_extract:
        result = CliRunner().invoke(
            precompute_cli, ["extract-query", str(sql_path), str(output_path)]
        )

    assert result.exit_code == 0, result.output
    mock_extract.assert_called_once_with(
        sql_path, output_path, "postgresql://user:pass@host/db", range=None
    )
