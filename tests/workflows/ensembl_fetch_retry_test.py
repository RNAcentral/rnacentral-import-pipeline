# -*- coding: utf-8 -*-

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENSEMBL = (ROOT / "workflows" / "databases" / "ensembl.nf").read_text()


def fetch_error_strategy() -> str:
    body = ENSEMBL.split("process fetch_species_data {", 1)[1]
    directives = body.split("input:", 1)[0]
    [line] = [
        l.strip()
        for l in directives.splitlines()
        if l.strip().startswith("errorStrategy")
    ]
    return line


# wget exits 4 on a network failure; under 'ignore' that silently drops the species.
def test_fetch_species_data_retries_network_failures():
    retried = re.search(r"exitStatus\s+in\s+\[([^\]]*)\]", fetch_error_strategy())
    assert retried and "4" in [s.strip() for s in retried.group(1).split(",")], (
        f"fetch_species_data does not retry wget exit 4: {fetch_error_strategy()!r}"
    )


def test_fetch_species_data_still_retries_server_errors():
    retried = re.search(r"exitStatus\s+in\s+\[([^\]]*)\]", fetch_error_strategy())
    assert retried and "8" in [s.strip() for s in retried.group(1).split(",")]


# A FULL release retires whatever is absent, so a species that never downloads must stop the run.
def test_fetch_species_data_never_drops_a_species_after_its_retries():
    assert "'ignore'" not in fetch_error_strategy()
    assert "'terminate'" in fetch_error_strategy()


def directives(path, process):
    body = (ROOT / path).read_text().split(f"process {process} {{", 1)[1]
    return [l.strip() for l in body.split("output:", 1)[0].splitlines() if l.strip()]


# Both read species.json first; the EBI FTP refusing connections for a minute killed a run.
def test_reading_species_json_retries_after_a_pause():
    for path, process in [
        ("workflows/databases/ensembl.nf", "find_urls"),
        ("workflows/metadata/ensembl.nf", "assemblies"),
    ]:
        lines = directives(path, process)
        [strategy] = [l for l in lines if l.startswith("errorStrategy")]
        assert "sleep(" in strategy and "'retry'" in strategy, (process, strategy)
        assert any(l.startswith("maxRetries") for l in lines), process
