# -*- coding: utf-8 -*-

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = (ROOT / "workflows" / "genome-mapping.nf").read_text()


def directives(process):
    body = WORKFLOW.split(f"process {process} {{", 1)[1]
    return body.split("input:", 1)[0]


# Both fetch from the new Ensembl FTP, which drops connections under load; 'ignore' alone skips the species.
@pytest.mark.parametrize("process", ["download_genome", "get_browser_coordinates"])
def test_downloads_retry_network_failures(process):
    strategy = re.search(r"errorStrategy\s*\{(.*)\}", directives(process))
    assert strategy, f"{process} has no retrying errorStrategy"
    retried = re.search(r"exitStatus\s+in\s+\[([^\]]*)\]", strategy.group(1))
    assert retried and "4" in [s.strip() for s in retried.group(1).split(",")]


def test_genomes_come_from_species_json_not_the_legacy_ftp():
    assert "url-for" not in WORKFLOW
    assert "update-assemblies" not in WORKFLOW
    assert "remote-files" in WORKFLOW
