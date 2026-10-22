# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
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

import json

import polars as pl
import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import cli
from rnacentral_pipeline.rnacentral.genes import convert
from rnacentral_pipeline.rnacentral.genes import lineage as L

GFF = """##gff-version 3
#!genome-version GRCh38
1\tRNAcentral\ttranscript\t100\t400\t.\t+\t.\tID=URS0000000001_9606@1/100-150%2C300-400:+;Name=URS0000000001_9606;type=lncRNA
1\tRNAcentral\texon\t300\t400\t.\t+\t.\tID=e2;Parent=URS0000000001_9606@1/100-150%2C300-400:+
1\tRNAcentral\texon\t100\t150\t.\t+\t.\tID=e1;Parent=URS0000000001_9606@1/100-150%2C300-400:+
1\tRNAcentral\ttranscript\t500\t600\t.\t-\t.\tID=t2;Name=URS0000000002;type=miRNA
1\tRNAcentral\ttranscript\t800\t900\t.\t+\t.\tID=t3;Name=URS0000000003;type=circRNA
1\tRNAcentral\tpredicted_gene\t90\t700\t.\t+\t.\tID=RNACENSG00000000001.1;Name=RNACENSG00000000001
"""


@pytest.fixture
def gff(tmp_path):
    path = tmp_path / "x.gff3"
    path.write_text(GFF)
    return path


def test_assembly_id_skips_the_regions_join_and_keeps_every_transcript(gff):
    df = convert.gff_to_polars(gff, 9606, assembly_id="GRCh38")
    assert df.height == 2
    assert set(df["assembly_id"]) == {"GRCh38"}


def test_region_name_is_taken_from_the_gff_id_and_rebuilt_otherwise(gff):
    df = convert.gff_to_polars(gff, 9606, assembly_id="GRCh38").sort("region_name")
    assert df["region_name"].to_list() == [
        "URS0000000001_9606@1/100-150,300-400:+",
        "URS0000000002_9606@1/500-600:-",
    ]


def test_predicted_genes_are_not_transcripts(gff):
    df = convert.gff_to_polars(gff, 9606, assembly_id="GRCh38")
    assert df["urs_taxid"].str.starts_with("URS").all()


def test_circrnas_are_left_out(gff):
    df = convert.gff_to_polars(gff, 9606, assembly_id="GRCh38")
    assert "URS0000000003_9606" not in set(df["urs_taxid"])


def test_get_assembly_has_no_trailing_newline(gff):
    assert convert.get_assembly(gff) == "GRCh38"


def test_needs_assembly_or_regions(gff):
    with pytest.raises(ValueError):
        convert.gff_to_polars(gff, 9606)


def transcripts():
    names = [
        "URS0000000001_9606@1/100-200:+",
        "URS0000000002_9606@1/150-250:+",
        "URS0000000003_9606@1/5000-5100:+",
    ]
    return pl.DataFrame(
        {
            "assembly_id": ["A1"] * 3,
            "chromosome": ["1"] * 3,
            "region_name": names,
            "region_start": [100, 150, 5000],
            "region_stop": [200, 250, 5100],
            "strand": [1, 1, 1],
        }
    )


def test_clusters_from_assignments_maps_strand_and_leaves_singletons_out():
    a = pl.DataFrame(
        {
            "region_name": transcripts()["region_name"][:2].to_list(),
            "gene_id": ["g1", "g1"],
        }
    )
    c = L.clusters_from_assignments(transcripts(), a)
    assert c.height == 2 and set(c["strand"]) == {"+"}
    assert set(c.columns) == set(L.CLUSTER_COLUMNS)
    L.build_clusters(c)  # accepted by the planner


def test_clusters_from_assignments_rejects_unknown_regions():
    a = pl.DataFrame({"region_name": ["nope@1/1-2:+"], "gene_id": ["g1"]})
    with pytest.raises(ValueError, match="not in the transcripts"):
        L.clusters_from_assignments(transcripts(), a)


def test_cli_chain_matches_calling_plan_directly(tmp_path):
    runner = CliRunner()
    t = transcripts()
    t.write_parquet(tmp_path / "t.pq")
    state = "-"
    direct = L.empty_state()
    for release, members in enumerate([[0, 1], [0, 1, 2]], start=1):
        names = t["region_name"][members].to_list()
        a = pl.DataFrame({"region_name": names, "gene_id": ["g"] * len(names)})
        a.write_parquet(tmp_path / f"a{release}.pq")
        out = tmp_path / f"r{release}"
        r = runner.invoke(
            cli,
            [
                "genes",
                "utils",
                "clusters",
                "--transcripts",
                str(tmp_path / "t.pq"),
                "--assignments",
                str(tmp_path / f"a{release}.pq"),
                "--output",
                str(tmp_path / f"c{release}.pq"),
            ],
        )
        assert r.exit_code == 0, r.output
        r = runner.invoke(
            cli,
            [
                "genes",
                "utils",
                "plan",
                "--state",
                state,
                "--clusters",
                str(tmp_path / f"c{release}.pq"),
                "--release",
                str(release),
                "--prefix",
                "TEST",
                "--out_dir",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        state = str(out / "state.parquet")
        direct = L.plan(
            direct,
            L.clusters_from_assignments(t, a),
            release,
            "TEST",
        )
        direct_state = L.as_state(direct.genes)
        assert pl.read_parquet(state).equals(direct_state)
        direct = direct_state
        assert json.loads((out / "stats.json").read_text())
    assert pl.read_parquet(tmp_path / "r2" / "events.parquet")[
        "event_type"
    ].to_list() == ["changed"]


def test_collect_gathers_a_chain_of_releases(tmp_path):
    runner = CliRunner()
    t = transcripts()
    t.write_parquet(tmp_path / "t.pq")
    state, dirs = "-", []
    for release, members in enumerate([[0, 1], [0, 1, 2]], start=1):
        names = t["region_name"][members].to_list()
        pl.DataFrame(
            {"region_name": names, "gene_id": ["g"] * len(names)}
        ).write_parquet(tmp_path / f"a{release}.pq")
        runner.invoke(
            cli,
            [
                "genes",
                "utils",
                "clusters",
                "--transcripts",
                str(tmp_path / "t.pq"),
                "--assignments",
                str(tmp_path / f"a{release}.pq"),
                "--output",
                str(tmp_path / f"c{release}.pq"),
            ],
        )
        out = tmp_path / f"r{release}"
        r = runner.invoke(
            cli,
            [
                "genes",
                "utils",
                "plan",
                "--state",
                state,
                "--clusters",
                str(tmp_path / f"c{release}.pq"),
                "--release",
                str(release),
                "--prefix",
                "TEST",
                "--out_dir",
                str(out),
            ],
        )
        assert r.exit_code == 0, r.output
        state = str(out / "state.parquet")
        dirs.append(str(out))
    final = tmp_path / "final"
    args = ["genes", "utils", "collect", "--output_dir", str(final)]
    for d in dirs:
        args += ["--plan_dir", d]
    r = runner.invoke(cli, args)
    assert r.exit_code == 0, r.output
    assert sorted(
        pl.read_parquet(final / "events.parquet")["release"].unique().to_list()
    ) == [1, 2]
    assert sorted(
        pl.read_parquet(final / "member_deltas.parquet")["release"].unique().to_list()
    ) == [1, 2]
    assert pl.read_parquet(final / "state.parquet").equals(
        pl.read_parquet(tmp_path / "r2" / "state.parquet")
    )
