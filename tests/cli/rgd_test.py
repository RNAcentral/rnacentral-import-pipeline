# -*- coding: utf-8 -*-

import obonet
from click.testing import CliRunner

# Stub the ontology load for the import only. Left in place, this leaked into
# every module imported after it in a full run: the real read_obo is called
# with ignore_obsolete, so a one argument lambda took out ~200 precompute tests
# with an arity error that only ever showed up in the whole suite.
_real_read_obo = obonet.read_obo
obonet.read_obo = lambda *args, **kwargs: {}
try:
    from rnacentral_pipeline.cli import rgd as rgd_cli
    from rnacentral_pipeline.databases import rgd as rgd_db
finally:
    obonet.read_obo = _real_read_obo


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
