from pathlib import Path
import sys
import types

import pytest

from rnacentral_pipeline.stopfree import scan
from rnacentral_pipeline.stopfree.data import StopfreeResult


FASTA_FIXTURE = Path(__file__).resolve().parents[2] / "data/attempted/rfam/raw.fasta"


@pytest.mark.stopfree
def test_read_fasta_records_normalizes_uracil():
    records = list(scan.read_fasta_records(FASTA_FIXTURE))

    assert records[0] == ("URS0000000019", "GAATACAACATTCTTTGACCCTGC")
    assert records[1][0] == "URS000000001B"


@pytest.mark.stopfree
def test_build_results_uses_probability_threshold():
    records = [("URS0000000001_9606", "ATGATGATG")]
    results = list(
        scan.build_results(
            records,
            [("URS0000000001_9606", 3)],
            [("URS0000000001_9606", 0.44)],
            [("URS0000000001_9606", 0.01)],
            max_probability=0.05,
        )
    )

    assert results == [
        StopfreeResult.build(
            "URS0000000001_9606",
            3,
            0.44,
            0.01,
            0.05,
        )
    ]


@pytest.mark.stopfree
def test_calculate_results_uses_stopfree_module(monkeypatch):
    fake_stopfree = types.SimpleNamespace(
        calculate_stop_free_runs_with_ids=lambda records: [(records[0][0], 7)],
        calculate_gc_content=lambda records: [(records[0][0], 0.5)],
        calculate_run_probability=lambda runs, gcs: [(runs[0][0], 0.001)],
    )
    monkeypatch.setitem(sys.modules, "stopfree", fake_stopfree)

    results = scan.calculate_results(
        [("URS0000A423A4/1-83", "ATGATGATGATG")],
        max_probability=0.05,
    )

    assert results == [
        StopfreeResult.build(
            "URS0000A423A4/1-83",
            7,
            0.5,
            0.001,
            0.05,
        )
    ]


@pytest.mark.stopfree
def test_run_writes_csv(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        scan,
        "calculate_results",
        lambda records, max_probability: [
            StopfreeResult.build(
                records[0][0], 4, 0.4, 0.001, max_probability
            )
        ],
    )
    output = tmp_path / "out"

    scan.run(FASTA_FIXTURE, output, 0.05)

    assert (output / "results.csv").read_text() == "URS0000000019,4,0.4,0.001,True\n"
