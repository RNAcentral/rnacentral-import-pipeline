from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COORDINATES = (ROOT / "workflows" / "export" / "ftp" / "coordinates.nf").read_text()


def directives(process):
    body = COORDINATES.split(f"process {process} {{", 1)[1].split("input:", 1)[0]
    return [line.strip() for line in body.splitlines() if line.strip()]


def test_igv_gff3_retries_with_more_resources():
    # A SLURM kill (exit 140) on a large assembly aborted release 28's export.
    lines = directives("generate_gff3_for_igv")
    memory = [line for line in lines if line.startswith("memory ")]
    time = [line for line in lines if line.startswith("time ")]

    assert "errorStrategy 'retry'" in lines
    assert len(memory) == 1 and "task.attempt" in memory[0]
    assert len(time) == 1 and "task.attempt" in time[0]
