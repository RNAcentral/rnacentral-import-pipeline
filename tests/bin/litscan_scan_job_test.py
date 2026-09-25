# -*- coding: utf-8 -*-

"""
Tests for the parts of ``bin/litscan-scan-job.py`` that need no database or
Slurm environment: XML body-section parsing and XML-shard selection.

Pinning two PR-review fixes:
- parse_body_to_dict used to assign the deduped name to the wrong variable
  (``section_name`` instead of ``key``), so a second "Results" section
  silently overwrote the first instead of getting a unique key.
- add_xml_shard used a strict ``>`` when picking a shard, so a PMCID sitting
  exactly on a shard's starting boundary fell back to the previous shard.
"""

import importlib.machinery
import importlib.util
from pathlib import Path

import pytest
from lxml import etree as ET

BIN = Path(__file__).resolve().parents[2] / "bin" / "litscan-scan-job.py"


def _load_script():
    spec = importlib.util.spec_from_loader(
        "litscan_scan_job",
        importlib.machinery.SourceFileLoader("litscan_scan_job", str(BIN)),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


litscan_scan_job = _load_script()


def test_parse_body_to_dict_dedupes_repeated_section_names():
    body = ET.fromstring(
        "<body>"
        "<sec><title>Results</title><p>first</p></sec>"
        "<sec><title>Results</title><p>second</p></sec>"
        "</body>"
    )

    sections = litscan_scan_job.parse_body_to_dict(body, ignore_tags=set())

    assert sections["results"] == "first"
    assert sections["results_2"] == "second"


@pytest.mark.parametrize(
    "pmcid,expected_shard",
    [
        ("PMC100", "shard_100.xml"),  # exactly on a shard boundary
        ("PMC150", "shard_100.xml"),  # between two boundaries
        ("PMC200", "shard_200.xml"),  # exactly on the higher boundary
        ("PMC50", "shard_100.xml"),  # below smallest shard -> falls back to it
    ],
)
def test_add_xml_shard_boundary_is_inclusive(pmcid, expected_shard):
    xml_lookup = {100: Path("shard_100.xml"), 200: Path("shard_200.xml")}

    result = litscan_scan_job.add_xml_shard(pmcid, xml_lookup)

    assert result == expected_shard
