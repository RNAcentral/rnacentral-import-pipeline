# -*- coding: utf-8 -*-

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = (ROOT / "workflows" / "genome-mapping.nf").read_text()


def directives(process):
    body = WORKFLOW.split(f"process {process} {{", 1)[1]
    return body.split("input:", 1)[0]


# The Ensembl FTP drops connections under load; /nfs has the same tree but only datamover nodes see it.
@pytest.mark.parametrize(
    "process", ["setup", "download_genome", "get_browser_coordinates"]
)
def test_ensembl_files_are_read_from_nfs_on_a_datamover(process):
    body = WORKFLOW.split(f"process {process} {{", 1)[1].split("\n}\n", 1)[0]
    assert "queue 'datamover'" in body
    assert "--bind /nfs:/nfs" in body
    assert "wget" not in body


def test_genomes_come_from_species_json_not_the_legacy_ftp():
    assert "url-for" not in WORKFLOW
    assert "update-assemblies" not in WORKFLOW
    assert "remote-files" in WORKFLOW
