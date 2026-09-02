# -*- coding: utf-8 -*-

"""
Postgres fixtures for the SQL tests.

These tests run the query files in files/ verbatim against a real server, so the
thing under test is the artifact the pipeline ships rather than a paraphrase of
it. Set RNACENTRAL_TEST_DB to point at a scratch database, otherwise an
ephemeral cluster is started with the local initdb/pg_ctl. If neither is
available the tests skip.

WARNING: the fixtures DROP and CREATE the tables named in schema.sql. Never
point RNACENTRAL_TEST_DB at a database you care about.
"""

import os
import re
import shutil
import socket
import subprocess
import tempfile
from pathlib import Path

import psycopg2
import pytest

SCHEMA = Path(__file__).parent / "schema.sql"


def _free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


@pytest.fixture(scope="session")
def pg_url():
    if "RNACENTRAL_TEST_DB" in os.environ:
        url = os.environ["RNACENTRAL_TEST_DB"]
        # Check while the target is still pristine: once the tests start they
        # insert their own rows, and a per-test check would flag those.
        conn = psycopg2.connect(url)
        try:
            _refuse_if_populated(conn)
        finally:
            conn.close()
        yield url
        return

    initdb = shutil.which("initdb")
    pg_ctl = shutil.which("pg_ctl")
    if not initdb or not pg_ctl:
        pytest.skip("No RNACENTRAL_TEST_DB and no local initdb/pg_ctl")

    root = Path(tempfile.mkdtemp(prefix="rnac-pg-"))
    data = root / "data"
    port = _free_port()
    try:
        subprocess.run(
            [initdb, "-D", str(data), "-U", "postgres", "--auth=trust", "--no-sync"],
            check=True,
            capture_output=True,
        )
        subprocess.run(
            [
                pg_ctl,
                "-D",
                str(data),
                "-l",
                str(root / "log"),
                "-o",
                f"-p {port} -c listen_addresses=127.0.0.1 -c fsync=off",
                "-w",
                "start",
            ],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as err:
        shutil.rmtree(root, ignore_errors=True)
        pytest.skip(f"Could not start postgres: {err.stderr.decode()[:400]}")

    try:
        yield f"postgresql://postgres@127.0.0.1:{port}/postgres"
    finally:
        subprocess.run(
            [pg_ctl, "-D", str(data), "-m", "immediate", "-w", "stop"],
            capture_output=True,
        )
        shutil.rmtree(root, ignore_errors=True)


def _managed_tables():
    """The tables schema.sql drops, read from the file so the two cannot drift."""
    return re.findall(r"DROP TABLE IF EXISTS (\w+)", SCHEMA.read_text())


def _refuse_if_populated(conn):
    """Abort rather than DROP tables in a database that holds real data.

    RNACENTRAL_TEST_DB is easy to point at a real instance by accident and the
    schema fixture drops every table it manages, so check first. Resolution goes
    through to_regclass without a schema qualifier so it follows the connection's
    search_path -- a production instance may not keep these in `public`, and
    looking only there would report "absent" for a table the DROP would still
    find and destroy.

    A scratch database is empty, so any row in any managed table is enough to
    conclude we have the wrong target.
    """
    with conn.cursor() as cur:
        for table in _managed_tables():
            cur.execute("SELECT to_regclass(%s)", (table,))
            resolved = cur.fetchone()[0]
            if resolved is None:
                continue
            cur.execute(f"SELECT 1 FROM {resolved} LIMIT 1")
            if cur.fetchone() is not None:
                raise RuntimeError(
                    f"Refusing to run: `{resolved}` already has rows. These "
                    "tests DROP and recreate every table in schema.sql. Point "
                    "RNACENTRAL_TEST_DB at an empty scratch database."
                )


@pytest.fixture
def db(pg_url):
    """A connection with the fixture schema freshly created and empty."""
    conn = psycopg2.connect(pg_url)
    conn.autocommit = True
    with conn.cursor() as cur:
        cur.execute(SCHEMA.read_text())
    yield conn
    conn.close()
