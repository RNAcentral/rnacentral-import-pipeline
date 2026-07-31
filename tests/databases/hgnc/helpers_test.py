import pytest
import requests

from rnacentral_pipeline.databases.hgnc import helpers

BAD_ID = "ENSG00000235934"


class FakeResponse:
    def __init__(self, status_code, records=None):
        self.status_code = status_code
        self._records = records or []

    def raise_for_status(self):
        if self.status_code >= 400:
            err = requests.exceptions.HTTPError(str(self.status_code))
            err.response = self
            raise err

    def json(self):
        return self._records


class FakeSession:
    """Ensembl 400s the whole POST if any single id in it is withdrawn."""

    def __init__(self, status_for_bad=400):
        self.status_for_bad = status_for_bad
        self.batches = []

    def post(self, _url, headers=None, json=None, timeout=None):
        ids = json["ids"]
        self.batches.append(list(ids))
        if BAD_ID in ids:
            return FakeResponse(self.status_for_bad)
        return FakeResponse(200, [{"query": i, "seq": "ACGT"} for i in ids])


def _use(monkeypatch, session):
    monkeypatch.setattr(helpers, "_ensembl_session", lambda: session)
    return session


def test_one_bad_id_does_not_lose_the_rest_of_the_batch(monkeypatch):
    session = _use(monkeypatch, FakeSession())

    records = helpers._fetch_batch(["a", BAD_ID, "c", "d"])

    # The withdrawn id is dropped; every other gene in the batch survives.
    assert sorted(r["query"] for r in records) == ["a", "c", "d"]
    # It bisected rather than retrying the whole batch forever.
    assert [BAD_ID] in session.batches


def test_bad_id_alone_is_skipped_not_raised(monkeypatch):
    _use(monkeypatch, FakeSession())
    assert helpers._fetch_batch([BAD_ID]) == []


@pytest.mark.parametrize("status", [429, 500, 503])
def test_transient_failure_raises_rather_than_dropping_genes(monkeypatch, status):
    # The session has already exhausted its retries by this point. Dropping the
    # batch would retire these xrefs by absence under an incremental load.
    _use(monkeypatch, FakeSession(status_for_bad=status))

    with pytest.raises(RuntimeError):
        helpers._fetch_batch(["a", BAD_ID, "c"])


def test_network_failure_raises(monkeypatch):
    class Exploding:
        def post(self, *_args, **_kwargs):
            raise requests.exceptions.ConnectionError("no route to host")

    _use(monkeypatch, Exploding())

    with pytest.raises(RuntimeError):
        helpers._fetch_batch(["a", "b"])
