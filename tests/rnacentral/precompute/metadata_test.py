# -*- coding: utf-8 -*-

from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli.precompute import cli as precompute_cli
from rnacentral_pipeline.rnacentral.precompute import metadata


def _write_parquet(path, rows):
    pl.DataFrame(rows).write_parquet(path)


def test_read_spine_reads_basic_columns_directly(tmp_path):
    raw_dir = tmp_path
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )

    spine = metadata.read_spine(raw_dir).collect()

    assert spine.to_dicts() == [
        {
            "id": 1,
            "urs_id": 1,
            "urs_taxid": "URS1_9606",
            "urs": "URS1",
            "taxid": 9606,
            "length": 100,
        }
    ]


def test_as_list_groups_multiple_raw_rows_into_a_list_of_structs(tmp_path):
    path = tmp_path / "coordinates.parquet"
    _write_parquet(path, [{"id": 1, "x": 10}, {"id": 1, "x": 11}, {"id": 2, "x": 20}])

    out = metadata.as_list(pl.scan_parquet(path), "coordinates").collect()

    by_id = {row["id"]: row["coordinates"] for row in out.to_dicts()}
    assert by_id[1] == [{"x": 10}, {"x": 11}]
    assert by_id[2] == [{"x": 20}]


def test_as_optional_struct_collapses_at_most_one_row_per_id(tmp_path):
    path = tmp_path / "r2dt-hits.parquet"
    _write_parquet(path, [{"id": 1, "model_id": 5540}])

    out = metadata.as_optional_struct(pl.scan_parquet(path), "r2dt_hits").collect()

    assert out.to_dicts() == [{"id": 1, "r2dt_hits": {"model_id": 5540}}]


def test_as_optional_scalar_extracts_and_renames_a_single_field(tmp_path):
    path = tmp_path / "stopfree.parquet"
    _write_parquet(path, [{"id": 1, "is_protein_coding": True}])

    out = metadata.as_optional_scalar(
        pl.scan_parquet(path), "is_protein_coding", "possible_orf_stopfree"
    ).collect()

    assert out.to_dicts() == [{"id": 1, "possible_orf_stopfree": True}]


def test_build_plan_left_joins_a_struct_attachment_and_leaves_missing_ids_null(
    tmp_path,
):
    raw_dir = tmp_path
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            },
            {
                "id": 2,
                "urs_id": 2,
                "urs_taxid": "URS2_9606",
                "urs": "URS2",
                "taxid": 9606,
                "length": 100,
            },
        ],
    )
    _write_parquet(raw_dir / "r2dt-hits.parquet", [{"id": 1, "model_id": 5540}])
    attachments = [
        metadata.Attachment(
            "r2dt-hits.parquet",
            lambda raw: metadata.as_optional_struct(raw, "r2dt_hits"),
        ),
    ]

    plan = metadata.build_plan(raw_dir, attachments).collect().sort("id")

    rows = plan.to_dicts()
    assert rows[0]["r2dt_hits"] == {"model_id": 5540}
    assert rows[1]["r2dt_hits"] is None


def test_build_plan_fills_missing_ids_with_empty_list_for_list_attachments(tmp_path):
    # Review Focus item: an id with ZERO raw rows in a List-cardinality
    # attachment must get [], not null - a naive left join gives null,
    # because the id never appears in the grouped attachment table at all.
    raw_dir = tmp_path
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            },
            {
                "id": 2,
                "urs_id": 2,
                "urs_taxid": "URS2_9606",
                "urs": "URS2",
                "taxid": 9606,
                "length": 100,
            },
        ],
    )
    _write_parquet(raw_dir / "coordinates.parquet", [{"id": 1, "x": 1}])
    attachments = [
        metadata.Attachment(
            "coordinates.parquet",
            lambda raw: metadata.as_list(raw, "coordinates"),
            fill_nulls={"coordinates": []},
        ),
    ]

    plan = metadata.build_plan(raw_dir, attachments).collect().sort("id")

    rows = plan.to_dicts()
    assert rows[0]["coordinates"] == [{"x": 1}]
    assert rows[1]["coordinates"] == []


def test_build_plan_fails_loudly_when_an_attachment_file_is_missing(tmp_path):
    raw_dir = tmp_path
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )
    attachments = [
        metadata.Attachment(
            "coordinates.parquet", lambda raw: metadata.as_list(raw, "coordinates")
        ),
    ]

    with pytest.raises(FileNotFoundError):
        metadata.build_plan(raw_dir, attachments).collect()


def test_a_ropeatable_row_for_a_zero_or_one_attachment_does_not_crash(tmp_path):
    # Review Focus item: Rust's group() would error loudly on >1 row for a
    # ZeroOrOne attachment; this design's .first() silently picks one.
    # Document the accepted behavior rather than leaving it unverified.
    path = tmp_path / "r2dt-hits.parquet"
    _write_parquet(path, [{"id": 1, "model_id": 5540}, {"id": 1, "model_id": 9999}])

    out = metadata.as_optional_struct(pl.scan_parquet(path), "r2dt_hits").collect()

    assert out.to_dicts()[0]["r2dt_hits"]["model_id"] in (5540, 9999)


def test_orf_join_sources_is_unique_non_null_set(tmp_path):
    path = tmp_path / "orfs.parquet"
    _write_parquet(
        path,
        [
            {"id": 1, "source": "cpat", "is_protein_coding": None},
            {"id": 1, "source": "cpat", "is_protein_coding": None},
            {"id": 1, "source": "cpc2", "is_protein_coding": None},
        ],
    )

    out = metadata.orf_join(pl.scan_parquet(path)).collect()

    assert set(out.to_dicts()[0]["orf_info"]["sources"]) == {"cpat", "cpc2"}


def test_orf_join_orf_info_is_none_when_every_source_is_null(tmp_path):
    # Review Focus item: rows exist for the id, but none carry a source -
    # must be None (matching Rust's sources.is_empty() -> None), not
    # {"sources": []}.
    path = tmp_path / "orfs.parquet"
    _write_parquet(path, [{"id": 1, "source": None, "is_protein_coding": None}])

    out = metadata.orf_join(pl.scan_parquet(path)).collect()

    assert out.to_dicts()[0]["orf_info"] is None


def test_orf_join_possible_orf_is_first_non_null_not_any_true(tmp_path):
    # Rust: orfs.iter().find_map(|orf| orf.is_protein_coding) - the FIRST
    # non-null value in raw row order, not an OR over all rows.
    path = tmp_path / "orfs.parquet"
    _write_parquet(
        path,
        [
            {"id": 1, "source": None, "is_protein_coding": None},
            {"id": 1, "source": None, "is_protein_coding": False},
            {"id": 1, "source": None, "is_protein_coding": True},
        ],
    )

    out = metadata.orf_join(pl.scan_parquet(path)).collect()

    assert out.to_dicts()[0]["possible_orf"] is False


def test_attachments_list_has_all_seven_non_spine_types():
    names = {a.filename for a in metadata.ATTACHMENTS}
    assert names == {
        "coordinates.parquet",
        "rfam-hits.parquet",
        "r2dt-hits.parquet",
        "previous.parquet",
        "orfs.parquet",
        "stopfree.parquet",
        "tcode.parquet",
    }


# One dummy row per attachment type, at a sentinel id (0) outside every
# test's actual range - standing in for "no real data relevant to this
# test" for attachments the test doesn't otherwise care about.
_DUMMY_ATTACHMENT_ROWS = {
    "coordinates": {"id": 0, "x": 0},
    "rfam-hits": {"id": 0, "x": 0},
    "r2dt-hits": {"id": 0, "model_id": 0},
    "previous": {"id": 0, "x": 0},
    "orfs": {"id": 0, "source": None, "is_protein_coding": None},
    "stopfree": {"id": 0, "is_protein_coding": None},
    "tcode": {"id": 0, "is_protein_coding": None},
}


def _write_unused_attachments(raw_dir, names):
    for name in names:
        _write_parquet(raw_dir / f"{name}.parquet", [_DUMMY_ATTACHMENT_ROWS[name]])


def test_read_ranges_parses_upi_min_ut_min_ut_max_rows(tmp_path):
    path = tmp_path / "urs_taxid.csv"
    path.write_text("1,1,2\n3,3,5\n")

    assert metadata.read_ranges(path) == [(1, 1, 2), (3, 3, 5)]


def test_build_writes_absolute_paths_to_the_manifest_even_with_a_relative_output_dir(
    tmp_path, monkeypatch
):
    # Final review C1: build() used to write output_dir/name straight into
    # the manifest verbatim - if output_dir itself was relative (as
    # precompute.nf's `rnac precompute metadata-build $urs_taxid . metadata_out`
    # passes it), the manifest held relative paths. Nextflow's `file()` in
    # the workflow body resolves a relative path against the launch
    # directory, not the task's work directory where the file actually
    # lives - every downstream process_range would fail to find it.
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )
    _write_unused_attachments(
        raw_dir,
        [
            "coordinates",
            "rfam-hits",
            "r2dt-hits",
            "previous",
            "orfs",
            "stopfree",
            "tcode",
        ],
    )
    ranges_path = tmp_path / "urs_taxid.csv"
    ranges_path.write_text("1,1,1\n")

    monkeypatch.chdir(tmp_path)
    metadata.build(ranges_path, raw_dir, Path("out"))  # relative output_dir

    manifest_rows = (tmp_path / "out" / "manifest.csv").read_text().strip().splitlines()
    written_path = manifest_rows[0].split(",", 1)[1]
    assert Path(written_path).is_absolute()
    assert Path(written_path).exists()


def test_build_writes_one_output_file_per_range_and_a_manifest(tmp_path):
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": i,
                "urs_id": i,
                "urs_taxid": f"URS{i}_9606",
                "urs": f"URS{i}",
                "taxid": 9606,
                "length": 100,
            }
            for i in range(1, 6)
        ],
    )
    _write_unused_attachments(
        raw_dir,
        [
            "coordinates",
            "rfam-hits",
            "r2dt-hits",
            "previous",
            "orfs",
            "stopfree",
            "tcode",
        ],
    )

    ranges_path = tmp_path / "urs_taxid.csv"
    ranges_path.write_text("10,1,2\n20,3,5\n")

    output_dir = tmp_path / "out"
    metadata.build(ranges_path, raw_dir, output_dir)

    manifest_rows = (output_dir / "manifest.csv").read_text().strip().splitlines()
    assert len(manifest_rows) == 2

    chunk1 = pl.read_parquet(output_dir / "metadata-10.parquet")
    chunk2 = pl.read_parquet(output_dir / "metadata-20.parquet")
    assert chunk1.height == 2  # ids 1,2
    assert chunk2.height == 3  # ids 3,4,5


def test_build_logs_a_warning_with_dropped_count_for_an_unmatched_attachment_id(
    tmp_path, caplog
):
    # Spec's Error Handling section: an attachment id absent from basic is
    # silently dropped by the left join (unlike Rust's loud dense-fill
    # error) - must still be visible via a logged count, not silent.
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )
    _write_parquet(
        raw_dir / "coordinates.parquet", [{"id": 1, "x": 1}, {"id": 99, "x": 2}]
    )
    _write_unused_attachments(
        raw_dir, ["rfam-hits", "r2dt-hits", "previous", "orfs", "stopfree", "tcode"]
    )

    ranges_path = tmp_path / "urs_taxid.csv"
    ranges_path.write_text("1,1,1\n")
    output_dir = tmp_path / "out"

    with caplog.at_level("WARNING"):
        metadata.build(ranges_path, raw_dir, output_dir)

    assert any(
        "coordinates.parquet" in r.message and "1" in r.message for r in caplog.records
    )


def test_build_produces_a_valid_empty_file_for_a_range_matching_no_ids(tmp_path):
    # Review Focus item: an empty chunk must not crash sink_ndjson.
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )
    _write_unused_attachments(
        raw_dir,
        [
            "coordinates",
            "rfam-hits",
            "r2dt-hits",
            "previous",
            "orfs",
            "stopfree",
            "tcode",
        ],
    )

    ranges_path = tmp_path / "urs_taxid.csv"
    ranges_path.write_text("99,500,600\n")  # no id in basic falls in this range

    output_dir = tmp_path / "out"
    metadata.build(ranges_path, raw_dir, output_dir)

    assert pl.read_parquet(output_dir / "metadata-99.parquet").height == 0


def test_build_logs_a_warning_for_duplicate_rows_in_a_zero_or_one_attachment(
    tmp_path, caplog
):
    # Final review I2: Rust's grouper.rs errors loudly on >1 raw row for a
    # ZeroOrOne attachment; this design's .first() silently keeps one
    # instead. Must at least be logged, not fully silent - hides real data
    # problems (e.g. an unexpected duplicate r2dt hit for one URS).
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )
    _write_parquet(
        raw_dir / "r2dt-hits.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "model_id": 1,
                "model_name": "a",
                "model_source": "a",
                "model_so_term": None,
                "sequence_coverage": None,
                "model_coverage": None,
                "sequence_basepairs": None,
                "model_basepairs": None,
            },
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "model_id": 2,
                "model_name": "b",
                "model_source": "b",
                "model_so_term": None,
                "sequence_coverage": None,
                "model_coverage": None,
                "sequence_basepairs": None,
                "model_basepairs": None,
            },
        ],
    )
    _write_unused_attachments(
        raw_dir, ["coordinates", "rfam-hits", "previous", "orfs", "stopfree", "tcode"]
    )

    ranges_path = tmp_path / "urs_taxid.csv"
    ranges_path.write_text("1,1,1\n")
    output_dir = tmp_path / "out"

    with caplog.at_level("WARNING"):
        metadata.build(ranges_path, raw_dir, output_dir)

    assert any(
        "r2dt-hits.parquet" in r.message and "1" in r.message for r in caplog.records
    )


def test_metadata_build_cli_wires_through_to_metadata_build(tmp_path):
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    _write_parquet(
        raw_dir / "basic.parquet",
        [
            {
                "id": 1,
                "urs_id": 1,
                "urs_taxid": "URS1_9606",
                "urs": "URS1",
                "taxid": 9606,
                "length": 100,
            }
        ],
    )
    _write_unused_attachments(
        raw_dir,
        [
            "coordinates",
            "rfam-hits",
            "r2dt-hits",
            "previous",
            "orfs",
            "stopfree",
            "tcode",
        ],
    )
    ranges_path = tmp_path / "urs_taxid.csv"
    ranges_path.write_text("1,1,1\n")
    output_dir = tmp_path / "out"

    result = CliRunner().invoke(
        precompute_cli,
        ["metadata-build", str(ranges_path), str(raw_dir), str(output_dir)],
    )

    assert result.exit_code == 0, result.output
    assert (output_dir / "metadata-1.parquet").exists()
