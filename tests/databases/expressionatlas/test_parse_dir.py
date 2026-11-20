# -*- coding: utf-8 -*-

"""
`rnac expressionatlas parse-dir`, as workflows/databases/expressionatlas.nf
runs it on one experiment directory. The configuration files are real
Expression Atlas ones (E-MTAB-513: baseline, E-GEOD-21860: differential); the
data files are empty placeholders because the parsers are stood in for here --
they have their own tests -- so this checks which files get dispatched where.
"""

import shutil
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

from rnacentral_pipeline.databases.expressionatlas import cli, configuration, parser

DATA = Path("data/expressionatlas")


def test_parse_config_reads_a_baseline_experiment():
    config = configuration.parse_config(DATA / "E-MTAB-513-configuration.xml")
    assert config.exp_type == "rnaseq_mrna_baseline"
    groups = config.analytics[0].assay_groups.assay_group
    assert len(groups) == 16
    assert (groups[0].id, groups[0].label, groups[0].assays) == (
        "g1",
        "leukocyte",
        ["ERR030875"],
    )


def test_parse_config_reads_a_differential_experiment():
    config = configuration.parse_config(DATA / "E-GEOD-21860-configuration.xml")
    assert config.exp_type == "rnaseq_mrna_differential"
    (contrast,) = config.analytics[0].contrasts.contrast
    assert (contrast.id, contrast.ref_group, contrast.test_group) == (
        "g1_g2",
        "g1",
        "g2",
    )


HIT = pl.DataFrame({"urs_taxid": ["URS0000000001_9606"], "experiment": ["E-TEST"]})


@pytest.fixture
def calls(monkeypatch):
    seen = []

    def stand_in(kind):
        def parse(data, sdrf, lookup):
            seen.append((kind, Path(data).name, Path(sdrf).name, Path(lookup).name))
            return HIT

        return parse

    monkeypatch.setattr(parser, "parse_baseline", stand_in("baseline"))
    monkeypatch.setattr(parser, "parse_differential", stand_in("differential"))
    return seen


def experiment(tmp_path, accession, *files):
    directory = tmp_path / accession
    directory.mkdir()
    shutil.copy(DATA / f"{accession}-configuration.xml", directory)
    for name in files:
        (directory / name).write_text("")
    return directory


def parse_dir(tmp_path, directory):
    lookup = tmp_path / "lookup.csv"
    lookup.write_text("")
    output = tmp_path / "hits.ndjson"
    result = CliRunner().invoke(
        cli.cli, ["parse-dir", str(directory), str(lookup), str(output)]
    )
    assert result.exit_code == 0, result.output
    text = output.read_text()
    return pl.read_ndjson(output) if text.strip() else pl.DataFrame()


def test_a_baseline_experiment_uses_the_gene_tpms(tmp_path, calls):
    directory = experiment(
        tmp_path,
        "E-MTAB-513",
        "E-MTAB-513-tpms.tsv",
        "E-MTAB-513-transcripts-tpms.tsv",  # transcript level: not wanted
        "E-MTAB-513.condensed-sdrf.tsv",
    )
    hits = parse_dir(tmp_path, directory)
    assert calls == [
        (
            "baseline",
            "E-MTAB-513-tpms.tsv",
            "E-MTAB-513.condensed-sdrf.tsv",
            "lookup.csv",
        )
    ]
    assert hits.to_dicts() == HIT.to_dicts()


def test_a_differential_experiment_uses_the_analytics(tmp_path, calls):
    directory = experiment(
        tmp_path,
        "E-GEOD-21860",
        "E-GEOD-21860-analytics.tsv",
        "E-GEOD-21860.condensed-sdrf.tsv",
    )
    parse_dir(tmp_path, directory)
    assert calls == [
        (
            "differential",
            "E-GEOD-21860-analytics.tsv",
            "E-GEOD-21860.condensed-sdrf.tsv",
            "lookup.csv",
        )
    ]


@pytest.mark.parametrize(
    "accession,files",
    [
        # no condensed SDRF
        ("E-MTAB-513", ["E-MTAB-513-tpms.tsv"]),
        # two analytics files: ambiguous
        (
            "E-GEOD-21860",
            [
                "E-GEOD-21860-A-analytics.tsv",
                "E-GEOD-21860-B-analytics.tsv",
                "E-GEOD-21860.condensed-sdrf.tsv",
            ],
        ),
    ],
)
def test_an_incomplete_experiment_writes_no_hits(tmp_path, calls, accession, files):
    hits = parse_dir(tmp_path, experiment(tmp_path, accession, *files))
    assert calls == []
    assert hits.height == 0


def test_a_directory_without_a_configuration_writes_no_hits(tmp_path, calls):
    directory = tmp_path / "E-MTAB-513"
    directory.mkdir()
    hits = parse_dir(tmp_path, directory)
    assert calls == []
    assert hits.height == 0
