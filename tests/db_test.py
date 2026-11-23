# -*- coding: utf-8 -*-

import psycopg2
import pytest

from rnacentral_pipeline import db


def test_connect_retries_transient_refusal(monkeypatch):
    """A refusal that clears on its own should not surface to the caller."""

    calls = []

    def flaky(config, **options):
        calls.append(config)
        if len(calls) < 3:
            raise psycopg2.OperationalError(
                "server sent an error response during SSL exchange"
            )
        return "conn"

    monkeypatch.setattr(psycopg2, "connect", flaky)
    monkeypatch.setattr("time.sleep", lambda _: None)

    assert db.connect("postgres://example") == "conn"
    assert len(calls) == 3


def test_connect_gives_up(monkeypatch):
    """A server that never comes back must still raise, not hang forever."""

    calls = []

    def always_down(config, **options):
        calls.append(config)
        raise psycopg2.OperationalError("sorry, too many clients already")

    monkeypatch.setattr(psycopg2, "connect", always_down)
    monkeypatch.setattr("time.sleep", lambda _: None)

    with pytest.raises(psycopg2.OperationalError):
        db.connect("postgres://example")
    assert len(calls) == 6


def test_run_query_keeps_commit_flag_out_of_query_params(monkeypatch):
    """commit_on_leave configures the connection; it must not reach execute()."""

    seen = {}

    class FakeCursor:
        def execute(self, query, params):
            seen["params"] = params

        def __iter__(self):
            return iter([(1,)])

    from contextlib import contextmanager

    @contextmanager
    def fake_cursor(url, commit_on_leave=True):
        seen["commit_on_leave"] = commit_on_leave
        yield FakeCursor()

    monkeypatch.setattr(db, "cursor", fake_cursor)

    rows = list(
        db.run_query("url", "SELECT 1 WHERE x LIKE 'a%'", commit_on_leave=False)
    )

    assert rows == [(1,)]
    assert seen == {"commit_on_leave": False, "params": {}}
