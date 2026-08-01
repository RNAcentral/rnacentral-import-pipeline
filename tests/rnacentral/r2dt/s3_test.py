# -*- coding: utf-8 -*-

"""Tests for rnacentral_pipeline.rnacentral.r2dt.s3 (SVG upload to S3)."""

import hashlib
import io

import pytest
from botocore.exceptions import ClientError, ReadTimeoutError

from rnacentral_pipeline.rnacentral.r2dt import s3


@pytest.mark.r2dt
def test_s3_key_matches_legacy_layout():
    # $env/URS/${urs:3:2}/${urs:5:2}/${urs:7:2}/${urs:9:2}/<urs>.svg.gz
    assert (
        s3.s3_key("URS0000F7F700", "prod")
        == "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz"
    )
    assert (
        s3.s3_key("URS0000ABCD12", "test")
        == "test/URS/00/00/AB/CD/URS0000ABCD12.svg.gz"
    )


@pytest.mark.r2dt
def test_urs_of_strips_suffix():
    assert s3.urs_of("/some/dir/URS0000F7F700.svg.gz") == "URS0000F7F700"
    with pytest.raises(ValueError):
        s3.urs_of("/some/dir/URS0000F7F700.svg")


@pytest.mark.r2dt
def test_s3_key_rejects_bad_urs():
    with pytest.raises(ValueError):
        s3.s3_key("nope", "prod")


@pytest.mark.r2dt
def test_looks_like_md5():
    assert s3._looks_like_md5("d41d8cd98f00b204e9800998ecf8427e")
    assert not s3._looks_like_md5("d41d8cd9-abc")  # multipart-style etag
    assert not s3._looks_like_md5("")


@pytest.mark.r2dt
def test_client_prefers_the_pipeline_credentials(monkeypatch):
    captured = {}
    monkeypatch.setattr(
        s3.boto3, "client", lambda name, **kwargs: captured.update(kwargs) or name
    )

    monkeypatch.delenv("S3_KEY", raising=False)
    monkeypatch.delenv("S3_SECRET", raising=False)
    s3.client()
    assert "aws_access_key_id" not in captured

    monkeypatch.setenv("S3_KEY", "key")
    monkeypatch.setenv("S3_SECRET", "secret")
    s3.client()
    assert captured["aws_access_key_id"] == "key"
    assert captured["aws_secret_access_key"] == "secret"
    assert captured["endpoint_url"] == s3.ENDPOINT


class _FakeS3:
    """Minimal stand-in that records puts and can be told to fail one key."""

    def __init__(self, fail_key=None, objects=None, fail_times=None, error=None):
        self.puts = {}
        self.fail_key = fail_key
        self.fail_times = fail_times
        self.error = error or RuntimeError("simulated S3 500")
        self.attempts = 0
        self.objects = objects if objects is not None else {}

    def put_object(self, Bucket, Key, Body, ContentMD5, ContentType):
        if Key == self.fail_key:
            self.attempts += 1
            if self.fail_times is None or self.attempts <= self.fail_times:
                raise self.error
        self.puts[Key] = Body
        return {"ETag": '"%s"' % hashlib.md5(Body).hexdigest()}

    def head_object(self, Bucket, Key):
        if Key not in self.objects:
            raise ClientError(
                {"Error": {"Code": "404", "Message": "Not Found"}}, "HeadObject"
            )
        body, etag = self.objects[Key]
        return {"ContentLength": len(body), "ETag": '"%s"' % etag}


def _write_list(tmp_path, urs_ids):
    for u in urs_ids:
        (tmp_path / f"{u}.svg.gz").write_bytes(b"gzipped-svg-for-" + u.encode())
    listing = tmp_path / "file-list"
    listing.write_text("\n".join(str(tmp_path / f"{u}.svg.gz") for u in urs_ids) + "\n")
    return listing


def _stored(tmp_path, urs):
    body = (tmp_path / f"{urs}.svg.gz").read_bytes()
    return body, hashlib.md5(body).hexdigest()


@pytest.mark.r2dt
def test_upload_puts_all_and_verifies_etag(tmp_path, monkeypatch):
    fake = _FakeS3()
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    n = s3.upload(str(listing), "prod", workers=2)

    assert n == 2
    assert "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz" in fake.puts
    assert "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz" in fake.puts


@pytest.mark.r2dt
def test_upload_raises_on_any_failure(tmp_path, monkeypatch):
    # The whole point: a failed upload must NOT be swallowed (the curl bug).
    fake = _FakeS3(fail_key="prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz")
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=2)


@pytest.mark.r2dt
def test_upload_retries_a_read_timeout(tmp_path, monkeypatch):
    # One transient timeout used to fail the whole publish_layout task, throwing
    # away ~10 minutes of publishing for 1 of 2660 objects.
    key = "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz"
    fake = _FakeS3(fail_key=key, fail_times=2, error=ReadTimeoutError(endpoint_url=""))
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    assert s3.upload(str(listing), "prod", workers=2) == 2
    assert key in fake.puts


@pytest.mark.r2dt
def test_upload_gives_up_on_a_persistent_timeout(tmp_path, monkeypatch):
    fake = _FakeS3(
        fail_key="prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz",
        error=ReadTimeoutError(endpoint_url=""),
    )
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=2)
    assert fake.attempts == s3.MAX_ATTEMPTS


class _LandsThenTimesOutS3(_FakeS3):
    """A store that writes the object but loses the response.

    What livingobjects did to URS0000112770: the PUT completed server-side and
    the read timed out waiting for the ack, so every retry re-uploaded and
    timed out identically.
    """

    def put_object(self, Bucket, Key, Body, ContentMD5, ContentType):
        if Key == self.fail_key:
            self.attempts += 1
            self.puts[Key] = Body
            self.objects[Key] = (Body, hashlib.md5(Body).hexdigest())
            raise ReadTimeoutError(endpoint_url="")
        return super().put_object(Bucket, Key, Body, ContentMD5, ContentType)


@pytest.mark.r2dt
def test_upload_accepts_an_object_that_landed_despite_a_timeout(tmp_path, monkeypatch):
    """
    The regression: one object timing out on all 5 attempts failed the whole
    task ("1/3124 uploads failed"), and the Nextflow retry then redid the entire
    ~10-minute publish loop only to time out on the same object again.
    """
    key = "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz"
    fake = _LandsThenTimesOutS3(fail_key=key)
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    assert s3.upload(str(listing), "prod", workers=2) == 2
    # Settled on the first failure, so no pointless re-PUTs of a stored object.
    assert fake.attempts == 1


@pytest.mark.r2dt
def test_upload_still_fails_when_the_object_is_not_there(tmp_path, monkeypatch):
    """A timeout with nothing stored is a real failure and must stay loud."""
    key = "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz"
    fake = _FakeS3(fail_key=key, error=ReadTimeoutError(endpoint_url=""))
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=2)
    assert fake.attempts == s3.MAX_ATTEMPTS


@pytest.mark.r2dt
def test_upload_rejects_a_partial_object_left_by_a_timeout(tmp_path, monkeypatch):
    """
    Truncated bytes on the server are not a settled upload -- the size check has
    to send it back through the retries rather than call it done.
    """
    key = "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz"
    fake = _FakeS3(
        fail_key=key,
        error=ReadTimeoutError(endpoint_url=""),
        objects={key: (b"trunc", hashlib.md5(b"trunc").hexdigest())},
    )
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=1)
    assert fake.attempts == s3.MAX_ATTEMPTS


@pytest.mark.r2dt
def test_a_settled_check_does_not_rescue_a_permanent_error(tmp_path, monkeypatch):
    """403 is not transient; a stored object must not paper over it."""
    key = "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz"
    body = b"gzipped-svg-for-URS0000F7F701"
    fake = _FakeS3(
        fail_key=key,
        error=ClientError({"ResponseMetadata": {"HTTPStatusCode": 403}}, "PutObject"),
        objects={key: (body, hashlib.md5(body).hexdigest())},
    )
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    listing = _write_list(tmp_path, ["URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=1)
    assert fake.attempts == 1


@pytest.mark.r2dt
@pytest.mark.parametrize(
    "error,transient",
    [
        (ReadTimeoutError(endpoint_url=""), True),
        (ClientError({"ResponseMetadata": {"HTTPStatusCode": 503}}, "PutObject"), True),
        (
            ClientError({"ResponseMetadata": {"HTTPStatusCode": 403}}, "PutObject"),
            False,
        ),
        (ValueError("bad etag"), False),
    ],
)
def test_only_transient_errors_are_retried(error, transient):
    assert s3._is_transient(error) is transient


@pytest.mark.r2dt
def test_verify_passes_when_everything_matches(tmp_path, monkeypatch):
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])
    objects = {
        s3.s3_key(urs, "prod"): _stored(tmp_path, urs)
        for urs in ("URS0000F7F700", "URS0000F7F701")
    }
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeS3(objects=objects))

    out = io.StringIO()
    assert s3.verify(str(listing), "prod", workers=2, out=out) == 0
    assert out.getvalue() == ""


@pytest.mark.r2dt
def test_verify_reports_missing_and_corrupt_objects(tmp_path, monkeypatch):
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701", "URS0000F7F702"])
    good = s3.s3_key("URS0000F7F700", "prod")
    corrupt = s3.s3_key("URS0000F7F701", "prod")
    body, _ = _stored(tmp_path, "URS0000F7F701")
    objects = {
        good: _stored(tmp_path, "URS0000F7F700"),
        # right size, wrong checksum
        corrupt: (body, hashlib.md5(b"something else entirely").hexdigest()),
        # URS0000F7F702 is absent
    }
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeS3(objects=objects))

    out = io.StringIO()
    assert s3.verify(str(listing), "prod", workers=2, out=out) == 2
    reported = {
        line.split("\t")[1]: line.split("\t")[0] for line in out.getvalue().splitlines()
    }
    assert reported[corrupt].startswith("CHECKSUM")
    assert reported[s3.s3_key("URS0000F7F702", "prod")] == "MISSING"
    assert good not in reported


@pytest.mark.r2dt
def test_verify_reports_a_size_mismatch(tmp_path, monkeypatch):
    listing = _write_list(tmp_path, ["URS0000F7F700"])
    key = s3.s3_key("URS0000F7F700", "prod")
    monkeypatch.setattr(
        s3, "client", lambda *a, **k: _FakeS3(objects={key: (b"short", "x" * 32)})
    )

    out = io.StringIO()
    assert s3.verify(str(listing), "prod", workers=1, out=out) == 1
    assert out.getvalue().startswith("SIZE")


class _FakeListingS3:
    """Serves list_objects_v2 pages out of a flat key list."""

    def __init__(self, keys):
        self.keys = keys

    def get_paginator(self, name):
        assert name == "list_objects_v2"
        return self

    def paginate(self, Bucket, Prefix, Delimiter=None):
        matching = [k for k in self.keys if k.startswith(Prefix)]
        if Delimiter is None:
            yield {"Contents": [{"Key": k} for k in matching]}
            return
        prefixes = set()
        for key in matching:
            rest = key[len(Prefix) :]
            if Delimiter in rest:
                prefixes.add(Prefix + rest.split(Delimiter, 1)[0] + Delimiter)
        yield {"CommonPrefixes": [{"Prefix": p} for p in sorted(prefixes)]}


@pytest.mark.r2dt
def test_list_svgs_writes_every_urs(monkeypatch):
    keys = [
        "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz",
        "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz",
        "prod/URS/00/00/AB/CD/URS0000ABCD12.svg.gz",
        "prod/URS/00/01/AB/CD/URS0001ABCD12.svg.gz",
        "prod/URS/00/00/F7/F7/URS0000F7F700.notes",  # not an SVG
        "test/URS/00/00/F7/F7/URS0000F7F709.svg.gz",  # a different env
    ]
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeListingS3(keys))

    out = io.StringIO()
    n = s3.list_svgs("prod", out, workers=2)

    assert n == 4
    assert set(out.getvalue().split()) == {
        "URS0000F7F700",
        "URS0000F7F701",
        "URS0000ABCD12",
        "URS0001ABCD12",
    }


@pytest.mark.r2dt
@pytest.mark.parametrize("depth", [0, 1, 2, 4])
def test_list_svgs_is_complete_at_any_depth(monkeypatch, depth):
    keys = [
        "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz",
        "prod/URS/00/01/AB/CD/URS0001ABCD12.svg.gz",
        "prod/URS/01/02/03/04/URS0102030405.svg.gz",
    ]
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeListingS3(keys))

    out = io.StringIO()
    assert s3.list_svgs("prod", out, workers=2, depth=depth) == 3
    assert len(out.getvalue().split()) == 3
