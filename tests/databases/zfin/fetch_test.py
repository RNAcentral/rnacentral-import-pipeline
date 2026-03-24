# -*- coding: utf-8 -*-

import io
import importlib
import json
import urllib.error as error

import pytest

zfin = importlib.import_module("rnacentral_pipeline.databases.zfin.fetch")


def test_retries_transient_url_errors(monkeypatch):
    calls = {"count": 0}
    payload = {
        "metaData": {"publications": ["PMID: 12345"]},
        "data": [{"soTermId": "SO:0000035"}],
    }

    def fake_urlopen(url, timeout=30):
        calls["count"] += 1
        if calls["count"] < 3:
            raise error.URLError("temporary failure in name resolution")
        return io.StringIO(json.dumps(payload))

    monkeypatch.setattr(zfin.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(zfin, "sleep", lambda _: None)

    data = zfin.fetch("https://example.org/zfin.json")

    assert calls["count"] == 3
    assert data["metaData"]["publications"] == ["PMID:12345"]
    assert data["data"] == [{"soTermId": "SO:0000035"}]


def test_does_not_retry_non_transient_http_errors(monkeypatch):
    calls = {"count": 0}

    def fake_urlopen(url, timeout=30):
        calls["count"] += 1
        raise error.HTTPError(url, 404, "Not Found", hdrs=None, fp=None)

    monkeypatch.setattr(zfin.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(zfin, "sleep", lambda _: None)

    with pytest.raises(error.HTTPError):
        zfin.fetch("https://example.org/zfin.json")

    assert calls["count"] == 1
