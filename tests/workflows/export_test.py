# -*- coding: utf-8 -*-

"""
Guards for the export entry point (export.nf) and the workflows it runs.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXPORT = (ROOT / "export.nf").read_text()
CONFIG = (ROOT / "config/export.config").read_text()
FTP = (ROOT / "workflows/export/ftp.nf").read_text()
EXPORT_NFS = sorted((ROOT / "workflows/export").glob("**/*.nf"))


def block(text: str, name: str) -> str:
    """The brace-delimited body of `name { ... }` in a config file."""
    start = re.search(rf"^\s*{name}\s*\{{", text, re.MULTILINE).end()
    depth = 1
    for i in range(start, len(text)):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        if depth == 0:
            return text[start:i]
    raise AssertionError(f"unterminated block {name}")


def test_sequence_search_is_enabled_by_default():
    # find_db_to_export has `when: params.export?.sequence_search?.run`, so an
    # unset flag silently skips the whole export.
    assert re.search(r"^\s*run\s*=\s*true", block(CONFIG, "sequence_search"), re.M)


def test_qc_waits_for_the_ftp_export():
    assert "emit:" in FTP.split("workflow ftp {", 1)[1].split("\nworkflow {", 1)[0]
    assert "ftp.out" in EXPORT


def test_active_sequences_dump_only_runs_when_ftp_or_huggingface_does():
    # An ungated call runs the slowest query even when nothing consumes it.
    assert re.search(
        r"if \(params\.export\.ftp\.run \|\| params\.export\.huggingface\.run\) \{\s*"
        r"active_json = active_sequences\(\)",
        EXPORT,
    )
    assert (
        "active_sequences"
        not in FTP.split("workflow ftp {", 1)[1].split("\nworkflow {", 1)[0]
    )


def test_no_export_retry_guard_depends_on_an_exit_status():
    # A SLURM kill reaches nextflow with no exit status, so these never fire.
    offenders = [
        f"{p.relative_to(ROOT)}:{n}"
        for p in EXPORT_NFS
        for n, line in enumerate(p.read_text().splitlines(), 1)
        if "errorStrategy" in line and "exitStatus" in line
    ]
    assert not offenders, offenders


def test_publish_folders_default_to_the_launch_directory():
    assert "baseDir" not in CONFIG
    for section, folder in [
        ("search", "text-search"),
        ("sequence_search", "sequence-search"),
    ]:
        assert f"{{launchDir}}/{folder}" in block(block(CONFIG, section), "publish")
    assert "${launchDir}/sankey" in block(CONFIG, "sankey")
    assert "${launchDir}/ftp" in block(CONFIG, "ftp")


def test_export_gates_each_section_on_its_run_flag():
    assert "params.only" not in EXPORT
    for key in ["search", "sequence_search", "sankey", "ftp", "huggingface"]:
        assert f"if (params.export.{key}.run)" in EXPORT


def test_ftp_config_has_no_run_flags_that_nothing_reads():
    # Only the section-level flag is read (by export.nf); sub-section ones were dead.
    assert len(re.findall(r"\brun\s*=", block(CONFIG, "ftp"))) == 1


def test_each_shared_query_is_run_from_one_place():
    # The per-database dump served both the FTP and sequence-search exports, once each.
    users = [p.name for p in EXPORT_NFS if "database-specific.sql" in p.read_text()]
    assert users == ["database-sequences.nf"], users


def test_igv_gff3_is_derived_from_the_gff3_not_rebuilt():
    body = (ROOT / "workflows/export/ftp/coordinates.nf").read_text()
    igv = body.split("process generate_gff3_for_igv {", 1)[1].split("\nprocess ", 1)[0]
    assert "as-gff3" not in igv
