# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import csv
import hashlib
import io

import pyarrow.parquet as pq
import pytest

from rnacentral_pipeline.databases.ensembl.metadata import compara


@pytest.fixture
def simple():
    with open("data/ensembl/compara/simple.fasta", "r") as raw:
        yield raw


@pytest.mark.ensembl
def test_can_produce_expected_data(simple):
    assert list(compara.data(simple)) == [
        (
            "707b0ba6ac8febddca12a3126ff8bcfd7fd2838411b4316235162c43d0854907",
            [
                "ENSMUST00000082483",
                "ENSRNOT00000078435",
                "ENSRNOT00000082313",
                "ENSRNOT00000082656",
                "ENSRNOT00000082946",
                "ENSRNOT00000083118",
                "ENSRNOT00000083605",
                "ENSRNOT00000086169",
                "ENSRNOT00000086419",
                "ENSRNOT00000090157",
            ],
        ),
        (
            "2be126e36b38bb326dcde5358cfd1d9b76b8b8e46ea2b9876dacbca8ef894943",
            [
                "ENSMUST00000179813",
                "ENSMUST00000180265",
                "ENSRNOT00000078860",
                "ENSRNOT00000079197",
                "ENSRNOT00000087954",
                "ENSRNOT00000088084",
                "ENSRNOT00000091107",
            ],
        ),
    ]


# write(): blocks split by "//" become (homology_group, ensembl_transcript)
# rows, the group being the sha256 of the block's sorted transcript ids.
# Production's ensembl_compara is empty today, so these check the rules directly.

ALIGNMENTS = """>ENST00000516494
ACGU-ACGU
>ENSMUST00000083556
ACGUAACGU
//
>ENSMUST00000175389
UUGC
>ENST00000362432
UUGC
>ENSDART00000116137
UUGA
//
"""


def group(*ids):
    digest = hashlib.sha256()
    for eid in sorted(ids):
        digest.update(eid.encode("utf-8"))
    return digest.hexdigest()


EXPECTED = {
    ("ENST00000516494", group("ENST00000516494", "ENSMUST00000083556")),
    ("ENSMUST00000083556", group("ENST00000516494", "ENSMUST00000083556")),
    (
        "ENSMUST00000175389",
        group("ENSMUST00000175389", "ENST00000362432", "ENSDART00000116137"),
    ),
    (
        "ENST00000362432",
        group("ENSMUST00000175389", "ENST00000362432", "ENSDART00000116137"),
    ),
    (
        "ENSDART00000116137",
        group("ENSMUST00000175389", "ENST00000362432", "ENSDART00000116137"),
    ),
}


def test_each_block_is_one_homology_group():
    groups = list(compara.data(io.StringIO(ALIGNMENTS)))
    assert len(groups) == 2
    assert groups[0] == (
        group("ENST00000516494", "ENSMUST00000083556"),
        ["ENSMUST00000083556", "ENST00000516494"],
    )


def test_the_group_does_not_depend_on_block_order():
    reordered = ALIGNMENTS.replace(
        ">ENST00000516494\nACGU-ACGU\n>ENSMUST00000083556\nACGUAACGU\n",
        ">ENSMUST00000083556\nACGUAACGU\n>ENST00000516494\nACGU-ACGU\n",
    )
    assert list(compara.data(io.StringIO(reordered))) == list(
        compara.data(io.StringIO(ALIGNMENTS))
    )


@pytest.mark.parametrize("fmt", ["csv", "parquet"])
def test_write_gives_one_row_per_transcript(tmp_path, monkeypatch, fmt):
    monkeypatch.setenv("RNAC_OUTPUT_FORMAT", fmt)
    out = tmp_path / f"compara.{fmt}"
    compara.write(io.StringIO(ALIGNMENTS), out)
    if fmt == "csv":
        with out.open() as raw:
            rows = {(eid, key) for key, eid in csv.reader(raw)}
    else:
        rows = {
            (r["ensembl_transcript"], r["homology_group"])
            for r in pq.read_table(out).to_pylist()
        }
    assert rows == EXPECTED


def test_write_to_a_handle_is_csv():
    out = io.StringIO()
    compara.write(io.StringIO(ALIGNMENTS), out)
    assert {
        (eid, key) for key, eid in csv.reader(io.StringIO(out.getvalue()))
    } == EXPECTED
