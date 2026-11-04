# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import ensembl


@pytest.mark.cli
@pytest.mark.db
@pytest.mark.network
def test_can_fetch_assemblies():
    runner = CliRunner()
    base = Path(os.curdir).absolute()
    with runner.isolated_filesystem():
        Path("ucsc.json").write_text(
            '{"ucscGenomes": {"hg38": {"description": "Dec. 2013 (GRCh38/hg38)"}}}'
        )
        assemblies = "loaded-assemblies.csv"
        cmd = [
            "assemblies",
            str(base / "files" / "import-data" / "ensembl" / "example-locations.json"),
            "ucsc.json",
            assemblies,
        ]
        result = runner.invoke(ensembl.cli, cmd)
        assert result.exit_code == 0
        with open(assemblies) as raw:
            assert len(raw.readlines()) >= 5000


def parse_plasmodium_with(gff_text):
    runner = CliRunner()
    base = Path(os.curdir).absolute()
    with runner.isolated_filesystem():
        Path("genes.gff3").write_text(gff_text)
        return runner.invoke(
            ensembl.cli,
            [
                "parse",
                "--family-file",
                str(base / "data" / "rfam" / "families.tsv"),
                str(base / "data" / "ensembl" / "Plasmodium_falciparum_3D7.ncrna.embl"),
                "genes.gff3",
                ".",
            ],
        )


# Plasmodium once parsed to nothing and the run carried on; a FULL release would retire it.
@pytest.mark.cli
def test_parse_fails_when_the_gff_has_ncrnas_but_none_are_parsed():
    gff = Path("data/ensembl/Plasmodium_falciparum_3D7.ncrna.gff3").read_text()
    result = parse_plasmodium_with(gff.replace("PF3D7_1425100.1", "PF3D7_9999999.1"))
    assert result.exit_code != 0


@pytest.mark.cli
def test_parse_allows_an_empty_result_when_the_gff_has_no_ncrnas():
    gff = Path("data/ensembl/Plasmodium_falciparum_3D7.ncrna.gff3").read_text()
    result = parse_plasmodium_with(gff.replace("\tncRNA_gene\t", "\tgene\t"))
    assert result.exit_code == 0


# Any other error was once caught as "no entries", leaving a partial parse and a passing task.
@pytest.mark.cli
def test_parse_fails_loudly_on_any_other_error(monkeypatch):
    def broken(*args, **kwargs):
        yield from ()
        raise ValueError("could not convert string to float: '1:6734..6802'")

    monkeypatch.setattr(ensembl.parser, "parse", broken)
    result = parse_plasmodium_with(
        "##gff-version 3\n#!genome-version GCA_000002765\n#!genome-build-accession GCA_000002765.3\n"
    )
    assert result.exit_code != 0
    assert "1:6734..6802" in str(result.exception)


# Some genesets (Thermococcus kodakarensis) list ncRNAs as genes only, with no transcript or exon lines.
@pytest.mark.cli
def test_parse_uses_the_gene_when_the_gff_has_no_transcript_for_an_ncrna():
    lines = (
        Path("data/ensembl/Plasmodium_falciparum_3D7.ncrna.gff3")
        .read_text()
        .splitlines()
    )
    gff = "\n".join(l for l in lines if "\tsnoRNA\t" not in l and "\texon\t" not in l)
    result = parse_plasmodium_with(gff)
    assert result.exit_code == 0, result.output
