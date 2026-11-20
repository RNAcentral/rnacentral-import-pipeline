import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_make_rust_only_moves_binaries_the_workspace_builds():
    moved = re.findall(r"target/release/(\S+) bin", (ROOT / "Makefile").read_text())
    members = re.findall(r'"utils/([^"]+)"', (ROOT / "Cargo.toml").read_text())
    assert moved
    assert set(moved) <= set(members)
