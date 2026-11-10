# -*- coding: utf-8 -*-

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENSEMBL = (ROOT / "workflows" / "databases" / "ensembl.nf").read_text()


# A FULL release retires whatever is absent, so a species that fails to fetch must stop the run.
def test_fetch_species_data_never_drops_a_species_when_it_fails():
    body = ENSEMBL.split("process fetch_species_data {", 1)[1]
    assert "ignore" not in body.split("input:", 1)[0]


def directives(path, process):
    body = (ROOT / path).read_text().split(f"process {process} {{", 1)[1]
    return [l.strip() for l in body.split("output:", 1)[0].splitlines() if l.strip()]


# The EBI FTP refusing connections for a minute killed a run; /nfs has the same tree without that.
def test_find_urls_reads_species_json_from_nfs_on_a_datamover():
    lines = directives("workflows/databases/ensembl.nf", "find_urls")
    assert "queue 'datamover'" in lines
    assert "--location '${params.databases.ensembl.nfs_path}/species.json'" in ENSEMBL


def test_assemblies_still_retries_reading_species_json_after_a_pause():
    lines = directives("workflows/metadata/ensembl.nf", "assemblies")
    [strategy] = [l for l in lines if l.startswith("errorStrategy")]
    assert "sleep(" in strategy and "'retry'" in strategy, strategy
    assert any(l.startswith("maxRetries") for l in lines)
