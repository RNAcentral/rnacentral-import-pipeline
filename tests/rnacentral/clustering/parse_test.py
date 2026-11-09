# -*- coding: utf-8 -*-

import re
from pathlib import Path

import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.clustering.parse import parse

SCHEMA = Path(__file__).parents[3] / "files" / "clustering" / "schema.sql"

MEMBER = {"parent_urs": ["A"], "parent_len": [100], "n_fragments": [1]}


def write_members(path, urs, identity, ratio=0.5):
    pl.DataFrame(
        {
            **MEMBER,
            "fragment_urs": [urs],
            "fragment_len": [[50]],
            "fragment_ratio": [[ratio]],
            "fragment_identity": [identity],
        },
        schema_overrides={"fragment_identity": pl.List(pl.Float32)},
    ).write_parquet(path)


@pytest.fixture
def src(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    pl.DataFrame(
        {
            "urs": ["A", "B", "C", "D"],
            "cluster_95": ["A", "A", "A", "D"],
            "cluster_90": ["A", "A", "A", "A"],
            "length": [100, 50, 99, 80],
            "so_types": ["SO:0000655", None, None, None],
            "is_representative_95": [True, False, False, True],
            "is_representative_90": [True, False, False, False],
        }
    ).write_parquet(src / "clusters.parquet")
    return src


def test_clusters_keep_only_the_loaded_columns(src, tmp_path):
    write_members(src / "fragments.parquet", ["B"], [0.97])
    parse(src, tmp_path / "out")
    clusters = pl.read_parquet(tmp_path / "out" / "clusters.parquet")
    assert clusters.columns == ["urs", "cluster_95", "cluster_90"]
    assert clusters.height == 4


def test_members_are_exploded_and_labelled_by_source_file(src, tmp_path):
    write_members(src / "fragments.parquet", ["B"], [0.97])
    write_members(src / "variants.parquet", ["C"], [None], ratio=0.99)
    parse(src, tmp_path / "out")

    absorbed = pl.read_parquet(tmp_path / "out" / "absorption.parquet").sort("urs")
    assert absorbed.rows() == [
        ("B", "A", "fragment", pytest.approx(0.97, abs=1e-6)),
        ("C", "A", "variant", None),
    ]


def test_works_with_only_one_member_file(src, tmp_path):
    write_members(src / "variants.parquet", ["C"], [0.9])
    parse(src, tmp_path / "out")
    absorbed = pl.read_parquet(tmp_path / "out" / "absorption.parquet")
    assert absorbed["absorbed_as"].to_list() == ["variant"]


def test_no_member_files_is_an_error(src, tmp_path):
    with pytest.raises(ValueError):
        parse(src, tmp_path / "out")


@pytest.mark.parametrize(
    "table,parquet",
    [
        ("load_sequence_cluster", "clusters.parquet"),
        ("load_sequence_absorption", "absorption.parquet"),
    ],
)
def test_parser_output_matches_the_staging_tables(src, tmp_path, table, parquet):
    # load-parquet inserts by column name, so a rename on either side breaks the load.
    write_members(src / "fragments.parquet", ["B"], [0.97])
    parse(src, tmp_path / "out")

    ddl = re.search(rf"{table} \((.*?)\);", SCHEMA.read_text(), re.S).group(1)
    columns = [line.split()[0] for line in ddl.strip().splitlines()]
    assert columns == pl.read_parquet(tmp_path / "out" / parquet).columns
