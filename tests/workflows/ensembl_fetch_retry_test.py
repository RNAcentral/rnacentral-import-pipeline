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
def test_species_json_is_read_from_nfs_on_a_datamover():
    for path, process in [
        ("workflows/databases/ensembl.nf", "find_urls"),
        ("workflows/metadata/ensembl.nf", "assemblies"),
    ]:
        assert "queue 'datamover'" in directives(path, process), process
        body = (ROOT / path).read_text().split(f"process {process} {{", 1)[1]
        body = body.split("\n}\n", 1)[0]
        assert "--location '${params.databases.ensembl.nfs_path}/species.json'" in body


# Over HTTPS each md5sum.txt took ~19 s from a compute node, so 11.6k of them took hours.
def test_checksums_are_read_from_nfs_on_a_datamover():
    assert "queue 'datamover'" in directives(
        "workflows/databases/ensembl.nf", "ensembl_source_diff"
    )
    body = ENSEMBL.split("process ensembl_source_diff {", 1)[1].split("\n}\n", 1)[0]
    assert "--nfs-path '${params.databases.ensembl.nfs_path}'" in body
