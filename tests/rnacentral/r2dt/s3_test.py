# -*- coding: utf-8 -*-

"""Tests for rnacentral_pipeline.rnacentral.r2dt.s3 (SVG upload to S3)."""

import gzip
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


@pytest.mark.r2dt
def test_botocore_does_not_retry_underneath_us(monkeypatch):
    """
    _with_retries gives every object MAX_ATTEMPTS goes. Letting botocore retry
    as well multiplied the budgets, so one wedged key took 5 x 5 x 120s to fail.
    """
    captured = {}
    monkeypatch.setattr(
        s3.boto3, "client", lambda name, **kwargs: captured.update(kwargs) or name
    )

    s3.client()

    assert captured["config"].retries["max_attempts"] == 1


class _FakeS3:
    """Minimal stand-in that records puts and can be told to fail one key."""

    def __init__(self, fail_key=None, objects=None, fail_times=None, error=None):
        self.puts = {}
        self.fail_key = fail_key
        self.fail_times = fail_times
        self.error = error or RuntimeError("simulated S3 500")
        self.attempts = 0
        self.objects = objects if objects is not None else {}
        self.deleted = []

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

    def get_object(self, Bucket, Key):
        if Key == self.fail_key:
            raise self.error
        if Key not in self.objects:
            raise ClientError(
                {"Error": {"Code": "404", "Message": "Not Found"}}, "GetObject"
            )
        body, _ = self.objects[Key]
        return {"Body": io.BytesIO(body)}

    def delete_objects(self, Bucket, Delete):
        keys = [o["Key"] for o in Delete["Objects"]]
        self.deleted.append(keys)
        errors = [
            {"Key": k, "Code": "AccessDenied"} for k in keys if k == self.fail_key
        ]
        return {"Deleted": [{"Key": k} for k in keys], "Errors": errors}


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
def test_upload_tolerates_failures_up_to_the_cap(tmp_path, monkeypatch):
    """
    A single key the store will not accept (every PUT and HEAD to it timing out)
    otherwise costs the whole batch, and after the retries the run.
    """
    key = "prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz"
    fake = _FakeS3(fail_key=key, error=ReadTimeoutError(endpoint_url=""))
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])
    failures = tmp_path / "upload-failures.txt"

    n = s3.upload(
        str(listing), "prod", workers=2, allow_failures=1, failure_list=str(failures)
    )

    assert n == 1  # the good object still landed
    assert "URS0000F7F701" in failures.read_text()
    assert "URS0000F7F700" not in failures.read_text()


@pytest.mark.r2dt
def test_upload_still_fails_above_the_cap(tmp_path, monkeypatch):
    """An outage must not slip through as tolerated breakage."""

    class _AllFail(_FakeS3):
        def put_object(self, Bucket, Key, Body, ContentMD5, ContentType):
            raise ReadTimeoutError(endpoint_url="")

    monkeypatch.setattr(s3, "client", lambda *a, **k: _AllFail())
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])
    failures = tmp_path / "upload-failures.txt"

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(
            str(listing),
            "prod",
            workers=2,
            allow_failures=1,
            failure_list=str(failures),
        )
    # Written before raising, so the record survives the failure.
    assert failures.exists()


@pytest.mark.r2dt
def test_upload_writes_no_failure_list_on_a_clean_run(tmp_path, monkeypatch):
    """The Nextflow output is `optional: true` and keys off this file existing."""
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeS3())
    listing = _write_list(tmp_path, ["URS0000F7F700"])
    failures = tmp_path / "upload-failures.txt"

    s3.upload(str(listing), "prod", workers=1, failure_list=str(failures))

    assert not failures.exists()


@pytest.mark.r2dt
def test_upload_defaults_to_tolerating_nothing(tmp_path, monkeypatch):
    """Strictness stays the default; tolerance is opt-in per run."""
    fake = _FakeS3(
        fail_key="prod/URS/00/00/F7/F7/URS0000F7F701.svg.gz",
        error=ReadTimeoutError(endpoint_url=""),
    )
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)
    listing = _write_list(tmp_path, ["URS0000F7F700", "URS0000F7F701"])

    with pytest.raises(RuntimeError, match="uploads failed"):
        s3.upload(str(listing), "prod", workers=2)


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


def _failures(tmp_path, *urs_ids):
    path = tmp_path / "upload-failures.txt"
    path.write_text("".join(f"{u}\tRead timeout\n" for u in urs_ids))
    return path


def _data_csv(tmp_path, *urs_ids):
    path = tmp_path / "data_1.csv"
    path.write_text("".join(f"{u},1,(((...))),0,3\n" for u in urs_ids))
    return path


def _attempted_parquet(tmp_path, *urs_ids):
    import pyarrow as pa
    import pyarrow.parquet as pq

    path = tmp_path / "attempted_1.parquet"
    pq.write_table(
        pa.table({"urs": list(urs_ids), "r2dt_version": ["2.0"] * len(urs_ids)}), path
    )
    return path


@pytest.mark.r2dt
def test_drop_failed_removes_the_urs_from_both_sides(tmp_path):
    """
    Both files must lose the row. Dropping only the hits would leave the URS in
    pipeline_tracking_traveler, and files/r2dt/setup.sql deletes tracked URS
    from the to-draw set -- so it would never be attempted again.
    """
    import pyarrow.parquet as pq

    failures = _failures(tmp_path, "URS0000112770")
    data = _data_csv(tmp_path, "URS0000112770", "URS0000F7F700")
    attempted = _attempted_parquet(tmp_path, "URS0000112770", "URS0000F7F700")

    dropped = s3.drop_failed(str(failures), [str(data), str(attempted)])

    assert dropped == 2
    assert "URS0000112770" not in data.read_text()
    assert "URS0000F7F700" in data.read_text()
    assert pq.read_table(attempted).column("urs").to_pylist() == ["URS0000F7F700"]


@pytest.mark.r2dt
def test_drop_failed_is_a_noop_without_failures(tmp_path):
    """The workflow passes the list unconditionally, placeholder and all."""
    data = _data_csv(tmp_path, "URS0000F7F700")
    before = data.read_text()

    empty = tmp_path / "empty.txt"
    empty.write_text("")

    assert s3.drop_failed(str(empty), [str(data)]) == 0
    assert s3.drop_failed(str(tmp_path / "missing.txt"), [str(data)]) == 0
    assert data.read_text() == before


@pytest.mark.r2dt
def test_failed_urs_skips_the_run_headers(tmp_path):
    """The published copy is appended to per run, one `#` header per run."""
    path = tmp_path / "upload-failures.txt"
    path.write_text(
        "# run_one 2026-08-01T10:00\n"
        "URS0000112770\tRead timeout\n"
        "# run_two 2026-08-02T10:00\n"
        "URS0000F7F701\tRead timeout\n"
    )

    assert s3.failed_urs(str(path)) == {"URS0000112770", "URS0000F7F701"}


@pytest.mark.r2dt
def test_drop_failed_caps_the_whole_run(tmp_path):
    """
    upload-s3 tolerates per batch, so 1 each across 1000 batches is 1000 lost
    SVGs. The aggregate list is where a run-wide limit can actually be enforced.
    """
    failures = _failures(tmp_path, "URS0000000001", "URS0000000002", "URS0000000003")
    data = _data_csv(tmp_path, "URS0000000001", "URS0000F7F700")

    with pytest.raises(RuntimeError, match="over the 2 allowed"):
        s3.drop_failed(str(failures), [str(data)], max_failures=2)

    # Raised before touching the files, so nothing is half-dropped.
    assert "URS0000000001" in data.read_text()

    assert s3.drop_failed(str(failures), [str(data)], max_failures=3) == 1


@pytest.mark.r2dt
def test_drop_failed_matches_the_urs_column_only(tmp_path):
    """A URS appearing in another field must not take the row with it."""
    failures = _failures(tmp_path, "URS0000112770")
    data = tmp_path / "data_1.csv"
    data.write_text("URS0000F7F700,1,URS0000112770,0,3\n")

    assert s3.drop_failed(str(failures), [str(data)]) == 0
    assert "URS0000F7F700" in data.read_text()


class _FakeListingS3:
    """Serves list_objects_v2 pages out of a flat key list.

    `page_size` splits the keys into continuation-token pages; `fail_after` makes
    the given number of pages raise once each, to exercise the retry.
    """

    def __init__(self, keys, page_size=None, fail_pages=()):
        self.keys = keys
        self.page_size = page_size
        self.fail_pages = set(fail_pages)
        self.calls = 0

    def list_objects_v2(self, Bucket, Prefix, Delimiter=None, ContinuationToken=None):
        self.calls += 1
        matching = sorted(k for k in self.keys if k.startswith(Prefix))

        if Delimiter is not None:
            # A delimited listing returns BOTH the sub-prefixes and the keys
            # sitting directly at this level, as the real API does. Returning
            # only CommonPrefixes made a whole class of short listing invisible
            # to these tests.
            prefixes = set()
            direct = []
            for key in matching:
                rest = key[len(Prefix) :]
                if Delimiter in rest:
                    prefixes.add(Prefix + rest.split(Delimiter, 1)[0] + Delimiter)
                else:
                    direct.append(key)
            return {
                "CommonPrefixes": [{"Prefix": p} for p in sorted(prefixes)],
                "Contents": [{"Key": k} for k in direct],
            }

        start = int(ContinuationToken or 0)
        if start in self.fail_pages:
            self.fail_pages.discard(start)  # transient: fails once, then serves
            raise ReadTimeoutError(endpoint_url="")

        size = self.page_size or len(matching) or 1
        chunk = matching[start : start + size]
        end = start + len(chunk)
        page = {"Contents": [{"Key": k} for k in chunk]}
        if end < len(matching):
            page["IsTruncated"] = True
            page["NextContinuationToken"] = str(end)
        return page


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


@pytest.mark.r2dt
def test_list_svgs_refuses_a_prefix_holding_both_keys_and_children(monkeypatch):
    """
    A level with sub-prefixes AND its own SVGs cannot be sharded: discovery
    replaces it with its children, so its own keys are listed by nobody, and the
    run still reports a clean total. _FakeListingS3 derives CommonPrefixes from
    the keys it holds, so it never produced this shape -- which is why the
    listing tests all passed while the inventory came back 68% short.
    """
    keys = [
        "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz",
        "prod/URS/00/00/URS0000000000.svg.gz",  # sits directly in 00/00/
    ]
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeListingS3(keys))

    # At depth 2 the mixed level is a leaf and gets scanned recursively, which
    # is correct. It is only lost once discovery expands past it.
    assert s3.list_svgs("prod", io.StringIO(), workers=2, depth=2) == 2

    with pytest.raises(RuntimeError, match="silently skip"):
        s3.list_svgs("prod", io.StringIO(), workers=2, depth=3)


class _DroppingS3(_FakeListingS3):
    """Returns an empty, untruncated page for `drop` the first `times` times.

    The livingobjects failure mode: under concurrency a leaf prefix answers with
    no Contents, no IsTruncated and no error, which is byte-for-byte what an
    empty prefix looks like.
    """

    def __init__(self, keys, drop, times=1):
        super().__init__(keys)
        self.drop = drop
        self.times = times

    def list_objects_v2(self, Bucket, Prefix, Delimiter=None, ContinuationToken=None):
        if Delimiter is None and Prefix == self.drop and self.times:
            self.times -= 1
            return {"Contents": []}
        return super().list_objects_v2(
            Bucket, Prefix, Delimiter=Delimiter, ContinuationToken=ContinuationToken
        )


@pytest.mark.r2dt
def test_list_svgs_retries_a_prefix_that_drops_its_keys(monkeypatch):
    keys = [
        "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz",
        "prod/URS/00/01/AB/CD/URS0001ABCD12.svg.gz",
    ]
    fake = _DroppingS3(keys, drop="prod/URS/00/01/", times=1)
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)

    out = io.StringIO()
    assert s3.list_svgs("prod", out, workers=2, depth=2) == 2

    listed = out.getvalue().split()
    assert len(listed) == len(set(listed)) == 2  # retried, not double-written


@pytest.mark.r2dt
def test_list_svgs_fails_when_a_prefix_never_yields_its_keys(monkeypatch):
    """
    Two consecutive production runs counted 10.8M and 13.6M of the same ~41.2M
    objects, both exiting clean, because a dropped response looks exactly like
    an empty prefix. A short inventory must be an error, not a number.
    """
    keys = [
        "prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz",
        "prod/URS/00/01/AB/CD/URS0001ABCD12.svg.gz",
    ]
    fake = _DroppingS3(keys, drop="prod/URS/00/01/", times=99)
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)

    with pytest.raises(RuntimeError, match="listed empty"):
        s3.list_svgs("prod", io.StringIO(), workers=2, depth=2)


@pytest.mark.r2dt
def test_list_svgs_allows_a_genuinely_empty_bucket(monkeypatch):
    # The root is the one prefix nobody discovered, so it is allowed to be empty.
    monkeypatch.setattr(s3, "client", lambda *a, **k: _FakeListingS3([]))

    assert s3.list_svgs("prod", io.StringIO(), workers=2) == 0


@pytest.mark.r2dt
def test_list_svgs_resumes_after_a_page_timeout(monkeypatch):
    """
    One read timeout on prod/URS/00/00/ used to kill an inventory run 2.6M
    objects in. The retry must resume from the continuation token, not restart
    the prefix: sync() COPYs this into a PRIMARY KEY table, so a duplicated URS
    would fail the load.
    """
    keys = [f"prod/URS/00/00/F7/F7/URS0000F7F7{i:02d}.svg.gz" for i in range(10)]
    fake = _FakeListingS3(keys, page_size=3, fail_pages=[3, 6])
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)

    out = io.StringIO()
    assert s3.list_svgs("prod", out, workers=2) == 10

    listed = out.getvalue().split()
    assert len(listed) == len(set(listed)) == 10


@pytest.mark.r2dt
def test_list_svgs_still_fails_on_a_persistent_timeout(monkeypatch):
    # A dead endpoint must not silently produce a short inventory. The gap lands
    # in sync()'s `missing` list (in the DB, seemingly not in S3), so the run
    # would re-render thousands of SVGs that are already in the bucket.
    keys = ["prod/URS/00/00/F7/F7/URS0000F7F700.svg.gz"]
    fake = _FakeListingS3(keys, page_size=1, fail_pages=[])
    fake.fail_pages = None  # every page raises

    def always_timeout(*a, **k):
        raise ReadTimeoutError(endpoint_url="")

    fake.list_objects_v2 = always_timeout
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)
    monkeypatch.setattr(s3.time, "sleep", lambda _: None)

    with pytest.raises(ReadTimeoutError):
        s3.list_svgs("prod", io.StringIO(), workers=2)


@pytest.mark.r2dt
def test_list_svgs_refuses_a_truncated_page_with_no_token(monkeypatch):
    """
    A store that says "more to come" but omits the continuation token must not
    end the prefix quietly: that reports a clean total for a short listing, and
    sync() would then read the whole gap as work to redo.
    """

    class _NoTokenS3:
        def list_objects_v2(self, Bucket, Prefix, Delimiter=None, **kwargs):
            if Delimiter is not None:
                return {"CommonPrefixes": []}
            return {
                "Contents": [{"Key": f"{Prefix}URS0000F7F700.svg.gz"}],
                "IsTruncated": True,  # ... and no NextContinuationToken
            }

    monkeypatch.setattr(s3, "client", lambda *a, **k: _NoTokenS3())

    with pytest.raises(RuntimeError, match="no NextContinuationToken"):
        s3.list_svgs("prod", io.StringIO(), workers=2)


class _FakeHeadS3:
    """HEADs succeed for the keys it holds, 404 otherwise."""

    def __init__(self, present_urs, env="prod"):
        self.keys = {s3.s3_key(u, env) for u in present_urs}

    def head_object(self, Bucket, Key):
        if Key not in self.keys:
            raise ClientError({"Error": {"Code": "404"}}, "HeadObject")
        return {"ContentLength": 1}


def _missing_file(tmp_path, urs_ids):
    path = tmp_path / "missing-svgs.txt"
    path.write_text("".join(f"{u}\n" for u in urs_ids))
    return str(path)


@pytest.mark.r2dt
def test_spot_check_passes_when_the_missing_set_is_real(tmp_path):
    urs = [f"URS0000F7F7{i:02d}" for i in range(50)]
    fake = _FakeHeadS3([])  # none of them are in the bucket: the list is honest

    assert s3.spot_check_missing(fake, "prod", _missing_file(tmp_path, urs)) == 0.0


@pytest.mark.r2dt
def test_spot_check_tolerates_uploads_that_raced_the_listing(tmp_path):
    # Objects written while the listing ran are legitimately in both places.
    urs = [f"URS0000F7F7{i:02d}" for i in range(100)]
    fake = _FakeHeadS3(urs[:3])

    rate = s3.spot_check_missing(fake, "prod", _missing_file(tmp_path, urs))
    assert rate == pytest.approx(0.03)


@pytest.mark.r2dt
def test_spot_check_refuses_a_missing_set_that_is_really_in_the_bucket(tmp_path):
    """
    The August 2026 shape: sync reported 29.3M sequences to redraw, 94% of which
    were in the bucket already. A short listing and an empty bucket produce the
    same numbers and the same exit code, so the only thing that separates them
    is asking the bucket.
    """
    urs = [f"URS0000F7F7{i:02d}" for i in range(100)]
    fake = _FakeHeadS3(urs[:94])

    with pytest.raises(RuntimeError, match="listing is short, not the bucket"):
        s3.spot_check_missing(fake, "prod", _missing_file(tmp_path, urs))


@pytest.mark.r2dt
def test_spot_check_samples_the_whole_file_not_just_the_head(tmp_path):
    # COPY emits the missing list in key order, so a head -n sample would only
    # ever see the low URS range.
    urs = [f"URS0000{i:06X}" for i in range(10000)]
    sample = s3._reservoir_sample(_missing_file(tmp_path, urs), 200)

    assert len(sample) == 200
    assert len(set(sample)) == 200
    assert max(urs.index(u) for u in sample) > 5000


def _bucket_with(urs_ids, env="prod"):
    """A fake bucket holding a gzipped SVG for each URS."""
    objects = {}
    for u in urs_ids:
        body = gzip.compress(f"<svg id='{u}'/>".encode())
        objects[s3.s3_key(u, env)] = (body, hashlib.md5(body).hexdigest())
    return objects


@pytest.mark.r2dt
def test_local_path_shards_like_the_bucket(tmp_path):
    assert s3.local_path(tmp_path, "URS0000F7F700") == (
        tmp_path / "URS/00/00/F7/F7/URS0000F7F700.svg"
    )


@pytest.mark.r2dt
def test_download_writes_ungzipped_files_and_a_manifest(tmp_path, monkeypatch):
    urs_ids = ["URS0000F7F700", "URS0001ABCD12"]
    fake = _FakeS3(objects=_bucket_with(urs_ids))
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)

    urs_list = tmp_path / "urs.txt"
    urs_list.write_text("\n".join(urs_ids) + "\n")
    dest = tmp_path / "out"

    assert s3.download(dest, "prod", urs_list=urs_list, workers=2) == 2

    svg = dest / "URS/00/00/F7/F7/URS0000F7F700.svg"
    assert svg.read_bytes() == b"<svg id='URS0000F7F700'/>"

    rows = {
        line.split("\t")[0]: line.split("\t")
        for line in (dest / "manifest.tsv").read_text().splitlines()
    }
    assert set(rows) == set(urs_ids)
    urs, relative, size, md5 = rows["URS0000F7F700"]
    assert dest / relative == svg
    assert int(size) == svg.stat().st_size
    assert md5 == hashlib.md5(svg.read_bytes()).hexdigest()


@pytest.mark.r2dt
def test_download_skips_what_is_already_on_disk(tmp_path, monkeypatch):
    urs_ids = ["URS0000F7F700", "URS0001ABCD12"]
    fake = _FakeS3(objects=_bucket_with(urs_ids))
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)

    urs_list = tmp_path / "urs.txt"
    urs_list.write_text("\n".join(urs_ids) + "\n")
    dest = tmp_path / "out"
    done = s3.local_path(dest, "URS0000F7F700")
    done.parent.mkdir(parents=True)
    done.write_bytes(b"already here")

    # Only the other URS is fetched, and the pre-existing file is left alone.
    assert s3.download(dest, "prod", urs_list=urs_list, workers=2) == 1
    assert done.read_bytes() == b"already here"
    assert (dest / "manifest.tsv").read_text().count("\n") == 1


@pytest.mark.r2dt
def test_download_raises_but_still_writes_the_rest(tmp_path, monkeypatch):
    urs_ids = ["URS0000F7F700", "URS0001ABCD12"]
    fake = _FakeS3(
        objects=_bucket_with(urs_ids),
        fail_key=s3.s3_key("URS0000F7F700", "prod"),
        error=ClientError(
            {
                "Error": {"Code": "AccessDenied"},
                "ResponseMetadata": {"HTTPStatusCode": 403},
            },
            "GetObject",
        ),
    )
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)

    urs_list = tmp_path / "urs.txt"
    urs_list.write_text("\n".join(urs_ids) + "\n")
    dest = tmp_path / "out"

    with pytest.raises(RuntimeError):
        s3.download(dest, "prod", urs_list=urs_list, workers=2)
    assert s3.local_path(dest, "URS0001ABCD12").exists()


@pytest.mark.r2dt
def test_delete_orphans_batches_at_the_api_limit(tmp_path, monkeypatch):
    fake = _FakeS3()
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)

    urs_ids = ["URS%010X" % n for n in range(2500)]
    orphans = tmp_path / "orphan-svgs.txt"
    orphans.write_text("\n".join(urs_ids) + "\n")

    assert s3.delete_orphans(orphans, "prod", workers=2) == 2500
    assert sorted(len(b) for b in fake.deleted) == [500, 1000, 1000]
    assert set(sum(fake.deleted, [])) == {s3.s3_key(u, "prod") for u in urs_ids}


@pytest.mark.r2dt
def test_delete_orphans_raises_on_a_reported_error(tmp_path, monkeypatch):
    fake = _FakeS3(fail_key=s3.s3_key("URS0000F7F700", "prod"))
    monkeypatch.setattr(s3, "client", lambda *a, **k: fake)

    orphans = tmp_path / "orphan-svgs.txt"
    orphans.write_text("URS0000F7F700\n")

    with pytest.raises(RuntimeError):
        s3.delete_orphans(orphans, "prod", workers=1)
