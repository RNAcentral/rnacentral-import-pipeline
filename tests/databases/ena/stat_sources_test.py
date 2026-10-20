import re
from pathlib import Path

ENA_NF = Path(__file__).resolve().parents[3] / "workflows" / "databases" / "ena.nf"


def render(script: str) -> str:
    """Groovy's escapes as Nextflow applies them: \\\\ -> \\, \\$ -> $, and an unknown
    escape such as \\( is dropped outright."""
    keep = {"\\": "\\", "$": "$", "t": "\t", "n": "\n"}
    return re.sub(r"\\(.)", lambda m: keep.get(m.group(1), ""), script)


def test_stat_sources_find_keeps_its_grouping_after_rendering():
    """
    The stat_sources signature used '\\( ... \\)' with single backslashes, which
    Nextflow rendered away. Without the grouping find's -printf binds only to the
    '*.tar' branch, so every .ncr.gz source hashed as empty and all 12,165 ENA
    sources came back with the same sha256('') signature.
    """
    match = re.search(r"signature=\\\$\((find -L .*?)\n", ENA_NF.read_text())
    assert match, "stat_sources find command not found"
    rendered = render(match.group(1))
    assert r"-type f \( -name '*.ncr.gz' -o -name '*.tar' \)" in rendered
