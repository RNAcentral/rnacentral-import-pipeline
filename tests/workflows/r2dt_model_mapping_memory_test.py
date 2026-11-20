from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R2DT = (ROOT / "workflows" / "r2dt.nf").read_text()


def test_fetch_model_mapping_is_not_capped_in_megabytes():
    # A 16 MB cap OOM-killed psql (exit 137) and aborted release 29's analyze run.
    body = R2DT.split("process fetch_model_mapping {", 1)[1].split("input:", 1)[0]
    memory = [line for line in body.splitlines() if line.strip().startswith("memory")]
    assert not any("MB" in line for line in memory)
