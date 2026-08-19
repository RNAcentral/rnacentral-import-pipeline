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
