# -*- coding: utf-8 -*-

import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.precompute import select_releases


def test_select_new_matches_the_rust_binarys_own_test_case():
    # Direct port of releases.rs's test_select_new: URS1/URS3 newer in
    # xref -> selected; URS2 equal -> skipped; URS6 only in xref (novel)
    # -> selected; URS7 only in known (gone from xref) -> skipped.
    xref_lf = pl.LazyFrame(
        {
            "urs_taxid": ["URS1_9606", "URS2_9606", "URS3_9606", "URS6_9606"],
            "xref_release": [604, 600, 605, 607],
        }
    )
    known_lf = pl.LazyFrame(
        {
            "urs_taxid": ["URS1_9606", "URS2_9606", "URS3_9606", "URS7_9606"],
            "known_release": [600, 600, 600, 600],
        }
    )

    selected = select_releases.select_new(xref_lf, known_lf)

    assert sorted(selected) == ["URS1_9606", "URS3_9606", "URS6_9606"]


def test_select_new_raises_when_a_known_release_is_newer_than_xref():
    # Direct port of releases.rs's test_select_new_known_newer - a genuine
    # data-integrity check, not a style choice.
    xref_lf = pl.LazyFrame({"urs_taxid": ["URS1_9606"], "xref_release": [600]})
    known_lf = pl.LazyFrame({"urs_taxid": ["URS1_9606"], "known_release": [700]})

    with pytest.raises(ValueError, match="URS1_9606"):
        select_releases.select_new(xref_lf, known_lf)


def test_select_new_selects_an_id_never_seen_in_known_at_all():
    # Review Focus: never-precomputed - known_release is null, not a
    # sentinel "release 0" that could accidentally compare wrong.
    xref_lf = pl.LazyFrame({"urs_taxid": ["URS9_9606"], "xref_release": [1]})
    known_lf = pl.LazyFrame(
        {"urs_taxid": [], "known_release": []},
        schema={"urs_taxid": pl.String, "known_release": pl.Int64},
    )

    assert select_releases.select_new(xref_lf, known_lf) == ["URS9_9606"]


def test_select_new_skips_an_id_no_longer_in_xref_at_all():
    # Review Focus: the (None, Some) case - present in known, entirely
    # gone from xref (not just deleted-but-present) - must be silently
    # skipped, not selected, not an error.
    xref_lf = pl.LazyFrame(
        {"urs_taxid": [], "xref_release": []},
        schema={"urs_taxid": pl.String, "xref_release": pl.Int64},
    )
    known_lf = pl.LazyFrame({"urs_taxid": ["URS9_9606"], "known_release": [5]})

    assert select_releases.select_new(xref_lf, known_lf) == []


def test_write_reads_both_parquet_files_and_writes_one_urs_taxid_per_line(tmp_path):
    xref_path = tmp_path / "xref.parquet"
    known_path = tmp_path / "known.parquet"
    pl.DataFrame(
        {"urs_taxid": ["URS1_9606", "URS2_9606"], "xref_release": [604, 600]}
    ).write_parquet(xref_path)
    pl.DataFrame({"urs_taxid": ["URS2_9606"], "known_release": [600]}).write_parquet(
        known_path
    )

    output_path = tmp_path / "urs.csv"
    select_releases.write(xref_path, known_path, output_path)

    assert output_path.read_text().splitlines() == ["URS1_9606"]


def test_select_outdated_cli_wires_through_to_write(tmp_path):
    from unittest.mock import patch

    from click.testing import CliRunner

    from rnacentral_pipeline.cli.precompute import cli as precompute_cli

    xref_path = tmp_path / "xref.parquet"
    known_path = tmp_path / "known.parquet"
    xref_path.touch()
    known_path.touch()
    output_path = tmp_path / "urs.csv"

    with patch("rnacentral_pipeline.cli.precompute.pre_select.write") as mock_write:
        result = CliRunner().invoke(
            precompute_cli,
            ["select-outdated", str(xref_path), str(known_path), str(output_path)],
        )

    assert result.exit_code == 0, result.output
    mock_write.assert_called_once_with(xref_path, known_path, output_path)
