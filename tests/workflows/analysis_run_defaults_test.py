import re
from pathlib import Path

import pytest

CONFIG = Path(__file__).resolve().parents[2] / "config"


@pytest.mark.parametrize(
    "name", ["cpat", "genome-mapping", "r2dt", "rfam", "stopfree", "tcode"]
)
def test_analyses_only_run_when_turned_on(name):
    text = (CONFIG / f"{name}.config").read_text()
    assert re.search(r"^\s*run\s*=\s*false\s*$", text, re.MULTILINE)
