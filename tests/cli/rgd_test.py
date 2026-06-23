# -*- coding: utf-8 -*-

from click.testing import CliRunner
import obonet

obonet.read_obo = lambda _: {}
from rnacentral_pipeline.cli import rgd as rgd_cli
from rnacentral_pipeline.databases import rgd as rgd_db


def test_can_parse_rgd_files(tmp_path, monkeypatch):
    monkeypatch.setattr(rgd_db.helpers.phy, "lineage", lambda _: "lineage")
    monkeypatch.setattr(rgd_db.helpers.phy, "common_name", lambda _: "Norway rat")
    monkeypatch.setattr(rgd_db.helpers.phy, "species", lambda _: "Rattus norvegicus")

    output = tmp_path / "out"
    output.mkdir()

    runner = CliRunner()
    result = runner.invoke(
        rgd_cli.cli,
        [
            "parse",
            "data/rgd/sequences.fa.gz",
            "data/rgd/rat_genes.txt",
            str(output),
        ],
    )

    assert result.exit_code == 0
    assert (output / "accessions.csv").exists()
