# -*- coding: utf-8 -*-

"""
rnc_import_tracker decides which databases the weekly run imports at all, so a run
that forgets rows silently re-imports everything -- or, worse once the table is
populated, skips what it should not.
"""

import pytest
from click.testing import CliRunner

from rnacentral_pipeline.cli import scan_imports

DATABASES = {"ENA": 1, "HGNC": 2, "ZFIN": 44}


class FakeCursor:
    """Records every statement, and answers the one SELECT the command makes."""

    def __init__(self, rows):
        self.rows = rows
        self.executed = []
        # execute_values reaches back through the cursor for the connection encoding.
        self.connection = FakeConnection(self)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def execute(self, statement, params=None):
        if isinstance(statement, bytes):
            statement = statement.decode()
        self.executed.append((statement, params))

    def mogrify(self, template, args):
        return (template % tuple(f"'{arg}'" for arg in args)).encode()

    def fetchall(self):
        return list(self.rows)

    def close(self):
        pass


class FakeConnection:
    encoding = "UTF8"

    def __init__(self, cursor):
        self._cursor = cursor
        self.committed = False

    def cursor(self, **kwargs):
        return self._cursor

    def commit(self):
        self.committed = True

    def close(self):
        pass


@pytest.fixture
def cursor(monkeypatch):
    cur = FakeCursor(sorted(DATABASES.items()))
    monkeypatch.setattr(
        scan_imports.psycopg2, "connect", lambda url: FakeConnection(cur)
    )
    return cur


def _run(tmp_path, cursor, lines):
    md5s = tmp_path / "latest_md5s.csv"
    md5s.write_text(lines)
    result = CliRunner().invoke(
        scan_imports.cli, ["update-tracker", str(md5s), "--db-url", "postgres://unused"]
    )
    assert result.exit_code == 0, result.output
    return " ".join(statement for statement, _ in cursor.executed)


def test_update_tracker_never_truncates_the_table(tmp_path, cursor):
    """
    The TRUNCATE sat inside the per-row insert loop, so each row wiped the one before
    it: the table held a single database, and every other one looked new for ever.
    """
    sql = _run(tmp_path, cursor, "ena,md5-ena\nhgnc,md5-hgnc\nzfin,md5-zfin\n")

    assert "TRUNCATE" not in sql


def test_update_tracker_records_every_database(tmp_path, cursor):
    sql = _run(tmp_path, cursor, "ena,md5-ena\nhgnc,md5-hgnc\nzfin,md5-zfin\n")

    for name, checksum in [
        ("ena", "md5-ena"),
        ("hgnc", "md5-hgnc"),
        ("zfin", "md5-zfin"),
    ]:
        assert f"'{name}'" in sql
        assert f"'{checksum}'" in sql


def test_update_tracker_replaces_only_the_databases_it_wrote(tmp_path, cursor):
    """A database absent from this run keeps whatever the tracker already said."""
    _run(tmp_path, cursor, "ena,md5-ena\n")

    deletes = [params for sql, params in cursor.executed if sql.startswith("DELETE")]
    assert deletes == [(["ena"],)]


def test_update_tracker_ignores_a_database_rnacentral_does_not_know(tmp_path, cursor):
    """The db_id is a foreign key, so an unknown name must not reach the insert."""
    sql = _run(tmp_path, cursor, "ena,md5-ena\nnot-a-database,md5-nope\n")

    assert "'ena'" in sql
    assert "not-a-database" not in sql
