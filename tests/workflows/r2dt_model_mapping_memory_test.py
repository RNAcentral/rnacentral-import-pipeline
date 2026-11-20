from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R2DT = (ROOT / "workflows" / "r2dt.nf").read_text()


def directives(process):
    body = R2DT.split(f"process {process} {{", 1)[1].split("input:", 1)[0]
    return [line.strip() for line in body.splitlines() if line.strip()]


def test_fetch_model_mapping_is_not_capped_in_megabytes():
    # A 16 MB cap OOM-killed psql (exit 137) and aborted release 29's analyze run.
    memory = [
        line for line in directives("fetch_model_mapping") if line.startswith("memory")
    ]
    assert not any("MB" in line for line in memory)


def test_extract_sequences_retries_with_more_memory():
    # prepare-sequences was OOM-killed at a flat 12GB in release 29's analyze run.
    lines = directives("extract_sequences")
    memory = [line for line in lines if line.startswith("memory")]
    retry = [line for line in lines if line.startswith("errorStrategy")]
    assert len(memory) == 1 and "task.attempt" in memory[0]
    assert len(retry) == 1 and "task.attempt" in retry[0]
