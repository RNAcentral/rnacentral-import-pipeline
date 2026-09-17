# -*- coding: utf-8 -*-

"""
The ENA snapshot hands us the occasional truncated ``.ncr.gz``. Feeding every
archive to one ``xargs -r zcat`` turned that into ``exit status 123`` for the
whole ``fetch_directory`` task -- a wgs batch is 125 directories and a two-day
rsync, all thrown away because eight files in ``jbu``/``jbv`` ended early.

Unpacking each archive on its own has to keep the readable ones, drop the
partial output of a bad one and finish with a zero exit.
"""

import gzip
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENA = (ROOT / "workflows" / "databases" / "ena.nf").read_text()


def fetch_directory_script() -> str:
    body = ENA.split("process fetch_directory {", 1)[1].split("\nprocess ", 1)[0]
    return body.split('"""', 1)[1].rsplit('"""', 1)[0]


def unpack_section() -> str:
    """The part of the script between the rsync and the chunking, as bash."""
    script = fetch_directory_script()
    start = script.index("find copied -type f -empty -delete")
    end = script.index("mkdir $name-chunks")
    return script[start:end].replace("\\$", "$").replace("${name}", "wgs")


def test_archives_are_not_unpacked_through_a_shared_xargs():
    script = fetch_directory_script()
    shared = re.findall(r"xargs[^\n|]*\b(?:zcat|gzip|tar)\b", script)
    assert not shared, (
        f"fetch_directory unpacks archives through xargs ({shared}): one truncated "
        "member makes xargs exit 123 and costs the whole batch"
    )


def gnu_zcat_env(tmp_path):
    """BSD zcat only looks for .Z, so stand GNU's behaviour in front of it."""
    stub = tmp_path / "bin"
    stub.mkdir()
    (stub / "zcat").write_text('#!/bin/sh\nexec gzip -dc "$@"\n')
    (stub / "zcat").chmod(stat.S_IRWXU)
    return dict(os.environ, PATH=f"{stub}{os.pathsep}{os.environ['PATH']}")


def test_a_truncated_archive_is_skipped_not_fatal(tmp_path):
    good = tmp_path / "copied" / "a"
    bad = tmp_path / "copied" / "b"
    good.mkdir(parents=True)
    bad.mkdir(parents=True)
    (good / "good.ncr.gz").write_bytes(gzip.compress(b"ID   GOOD\n//\n"))
    (bad / "good2.ncr.gz").write_bytes(gzip.compress(b"ID   ALSO_GOOD\n//\n"))
    (bad / "truncated.ncr.gz").write_bytes(
        gzip.compress(b"ID   TRUNCATED\n" + b"x" * 4096 + b"\n//\n")[:32]
    )
    (bad / "empty.ncr.gz").write_bytes(b"")

    run = subprocess.run(
        ["bash", "-ue", "-c", unpack_section()],
        cwd=tmp_path,
        env=gnu_zcat_env(tmp_path),
        capture_output=True,
        text=True,
    )

    assert run.returncode == 0, run.stderr
    assert (tmp_path / "wgs.ncr").read_text() == "ID   GOOD\n//\nID   ALSO_GOOD\n//\n"
    assert "truncated.ncr.gz" in run.stderr


def test_a_truncated_tar_is_skipped_not_fatal(tmp_path):
    if not shutil.which("tar"):
        return

    copied = tmp_path / "copied" / "a"
    copied.mkdir(parents=True)
    member = copied / "member.ncr.gz"
    member.write_bytes(gzip.compress(b"ID   FROM_TAR\n//\n"))
    subprocess.run(
        ["tar", "-cf", str(copied / "good.tar"), "-C", str(copied), "member.ncr.gz"],
        check=True,
    )
    (copied / "broken.tar").write_bytes((copied / "good.tar").read_bytes()[:512])
    member.unlink()

    run = subprocess.run(
        ["bash", "-ue", "-c", unpack_section()],
        cwd=tmp_path,
        env=gnu_zcat_env(tmp_path),
        capture_output=True,
        text=True,
    )

    assert run.returncode == 0, run.stderr
    assert "broken.tar" in run.stderr
    assert (tmp_path / "wgs.ncr").read_text() == "ID   FROM_TAR\n//\n"
