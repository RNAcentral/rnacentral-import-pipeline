# -*- coding: utf-8 -*-

"""Tests for rnacentral_pipeline.rnacentral.r2dt.s3 (SVG upload to S3)."""

import io

import pytest

from rnacentral_pipeline.rnacentral.r2dt import s3


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


def test_urs_of_strips_suffix():
    assert s3.urs_of("/some/dir/URS0000F7F700.svg.gz") == "URS0000F7F700"
    with pytest.raises(ValueError):
        s3.urs_of("/some/dir/URS0000F7F700.svg")


def test_s3_key_rejects_bad_urs():
    with pytest.raises(ValueError):
        s3.s3_key("nope", "prod")


def test_looks_like_md5():
    assert s3._looks_like_md5("d41d8cd98f00b204e9800998ecf8427e")
    assert not s3._looks_like_md5("d41d8cd9-abc")  # multipart-style etag
    assert not s3._looks_like_md5("")


class _FakeS3:
    """Minimal stand-in that records puts and can be told to fail one key."""

    def __init__(self, fail_key=None):
        self.puts = {}
        self.fail_key = fail_key

    def put_object(self, Bucket, Key, Body, ContentMD5, ContentType):
        if Key == self.fail_key:
            raise RuntimeError("simulated S3 500")
        import hashlib

        self.puts[Key] = Body
        return {"ETag": '"%s"' % hashlib.md5(Body).hexdigest()}


def _write_list(tmp_path, urs_ids):
    for u in urs_ids:
        (tmp_path / f"{u}.svg.gz").write_bytes(b"gzipped-svg-for-" + u.encode())
    listing = tmp_path / "file-list"
    listing.write_text("\n".join(str(tmp_path / f"{u}.svg.gz") for u in urs_ids) + "\n")
    return listing


def test_upload_puts_all_and_verifies_etag(tmp_path, monkeypatch):
    fake = _FakeS3()
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    n = s3.upload(str(listing), "prod", workers=2)

    assert n == 2
    assert "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz" in fake.puts
    assert "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz" in fake.puts


def test_upload_raises_on_any_failure(tmp_path, monkeypatch):
    # The whole point: a failed upload must NOT be swallowed (the curl bug).
    fake = _FakeS3(fail_key="prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz")
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=2)
