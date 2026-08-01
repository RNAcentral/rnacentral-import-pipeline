# -*- coding: utf-8 -*-

"""
S3 handling for r2dt SVGs.

Replaces the old bin/update-svg.sh, which used `curl` without --fail: curl
returns exit 0 on an HTTP 4xx/5xx, so upload errors were swallowed and the
Nextflow process reported success while nothing landed. That silently diverged
the S3 SVG set from the DB should_show set by >1M objects.

Everything here fails loudly instead: boto3 raises on any error, every PUT
carries a Content-MD5 so the server rejects a corrupt body, and upload() re-raises
if any single object failed so the pipeline step goes red.

Objects live at (matching the historic layout in update-svg.sh):
    {env}/URS/{urs[3:5]}/{urs[5:7]}/{urs[7:9]}/{urs[9:11]}/{urs}.svg.gz
"""

import base64
import hashlib
import logging
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3
from botocore.config import Config
from botocore.exceptions import BotoCoreError, ClientError

LOGGER = logging.getLogger(__name__)

ENDPOINT = "https://livingobjects.ebi.ac.uk"
BUCKET = "ebi-rnacentral"
CONTENT_TYPE = "application/octet-stream"  # matches the existing objects
SUFFIX = ".svg.gz"

# botocore's retries share one adaptive quota per client, so a burst of slow
# responses can exhaust it; retry each object here too, on its own backoff.
MAX_ATTEMPTS = 5
BACKOFF = 2.0


def client(endpoint=ENDPOINT, pool_connections=64):
    """A path-style boto3 S3 client for the Weka endpoint.

    Honours S3_KEY/S3_SECRET (as update-svg.sh did) first, else falls back to
    boto3's normal chain (~/.aws/credentials, AWS_PROFILE, instance role, ...).
    """
    kwargs = dict(
        endpoint_url=endpoint,
        config=Config(
            s3={"addressing_style": "path"},
            max_pool_connections=pool_connections,
            connect_timeout=30,
            read_timeout=120,
            retries={"max_attempts": 5, "mode": "adaptive"},
        ),
    )
    if os.environ.get("S3_KEY") and os.environ.get("S3_SECRET"):
        kwargs["aws_access_key_id"] = os.environ["S3_KEY"]
        kwargs["aws_secret_access_key"] = os.environ["S3_SECRET"]
    return boto3.client("s3", **kwargs)


def urs_of(path):
    """The bare URS for a prepared file, e.g. .../URS0000F7F700.svg.gz -> URS0000F7F700."""
    name = Path(path).name
    if not name.endswith(SUFFIX):
        raise ValueError(f"Not a prepared SVG file: {path}")
    return name[: -len(SUFFIX)]


def s3_key(urs, env):
    """The S3 key for a URS, byte-for-byte compatible with update-svg.sh.

    update-svg.sh built: $env/${urs:0:3}/${urs:3:2}/${urs:5:2}/${urs:7:2}/${urs:9:2}/<file>
    where ${urs:0:3} is always "URS".
    """
    if len(urs) < 11 or not urs.startswith("URS"):
        raise ValueError(f"Unexpected URS id: {urs!r}")
    return f"{env}/URS/{urs[3:5]}/{urs[5:7]}/{urs[7:9]}/{urs[9:11]}/{urs}{SUFFIX}"


def _md5(data: bytes):
    """Return (hex, base64) MD5 of the bytes. hex compares to the S3 ETag; base64
    goes in the Content-MD5 header so the server verifies on receipt."""
    digest = hashlib.md5(data).digest()
    return digest.hex(), base64.b64encode(digest).decode("ascii")


def _read_file_list(file_list):
    with open(file_list) as fh:
        return [line.strip() for line in fh if line.strip()]


def _is_transient(err):
    """True for errors worth another go: timeouts, dropped connections, 5xx, 429."""
    if isinstance(err, ClientError):
        status = err.response.get("ResponseMetadata", {}).get("HTTPStatusCode", 0)
        return status >= 500 or status == 429
    return isinstance(err, BotoCoreError)


def _with_retries(fn, key, attempts=MAX_ATTEMPTS, sleep=time.sleep):
    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except Exception as err:  # noqa: BLE001 - classified by _is_transient
            if attempt == attempts or not _is_transient(err):
                raise
            LOGGER.warning(
                "Retrying %s after %s (attempt %d/%d)", key, err, attempt, attempts
            )
            sleep(BACKOFF * 2 ** (attempt - 1))


def upload(file_list, env, endpoint=ENDPOINT, bucket=BUCKET, workers=32):
    """Upload every gzipped SVG named in `file_list` to S3, in parallel.

    Each PUT sends Content-MD5 (server-side integrity) and we assert the returned
    ETag matches (client-side). Raises RuntimeError if ANY object failed, so the
    caller/pipeline sees the failure rather than a silent partial upload.
    """
    paths = _read_file_list(file_list)
    s3 = client(endpoint, pool_connections=workers + 8)
    lock = threading.Lock()
    failures = []
    done = {"n": 0}

    def put_one(path):
        data = Path(path).read_bytes()
        hex_md5, b64_md5 = _md5(data)
        key = s3_key(urs_of(path), env)
        resp = _with_retries(
            lambda: s3.put_object(
                Bucket=bucket,
                Key=key,
                Body=data,
                ContentMD5=b64_md5,
                ContentType=CONTENT_TYPE,
            ),
            key,
        )
        etag = resp.get("ETag", "").strip('"')
        # Single-part PUT => ETag is the hex MD5. If the store doesn't do that we
        # still had server-side Content-MD5 verification, so only warn.
        if etag and _looks_like_md5(etag) and etag != hex_md5:
            raise ValueError(f"ETag {etag} != md5 {hex_md5} for {key}")
        return key

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(put_one, p): p for p in paths}
        for fut in as_completed(futures):
            path = futures[fut]
            try:
                fut.result()
            except Exception as err:  # noqa: BLE001 - collect, report all at end
                with lock:
                    failures.append((path, err))
                LOGGER.error("Failed to upload %s: %s", path, err)
            else:
                with lock:
                    done["n"] += 1

    LOGGER.info("Uploaded %d/%d objects to %s", done["n"], len(paths), bucket)
    if failures:
        raise RuntimeError(
            f"{len(failures)}/{len(paths)} uploads failed; first: "
            f"{failures[0][0]}: {failures[0][1]}"
        )
    return done["n"]


def _looks_like_md5(etag):
    return len(etag) == 32 and all(c in "0123456789abcdef" for c in etag.lower())


def verify(
    file_list, env, endpoint=ENDPOINT, bucket=BUCKET, workers=32, out=sys.stdout
):
    """Check every file in `file_list` against the bucket by checksum.

    For each local file, HEAD the corresponding key and compare the ETag to the
    local MD5 (and size). Prints one line per problem to `out` and returns the
    number of mismatches (0 == all good), so a caller can exit non-zero.
    """
    from botocore.exceptions import ClientError

    paths = _read_file_list(file_list)
    s3 = client(endpoint, pool_connections=workers + 8)
    lock = threading.Lock()
    problems = {"n": 0}

    def check(path):
        data = Path(path).read_bytes()
        hex_md5, _ = _md5(data)
        key = s3_key(urs_of(path), env)
        try:
            head = s3.head_object(Bucket=bucket, Key=key)
        except ClientError as err:
            code = err.response.get("Error", {}).get("Code")
            reason = "MISSING" if code in ("404", "NoSuchKey") else f"ERROR {code}"
            return path, key, reason
        etag = head.get("ETag", "").strip('"')
        if head.get("ContentLength") != len(data):
            return path, key, f"SIZE {head.get('ContentLength')} != {len(data)}"
        if _looks_like_md5(etag):
            if etag != hex_md5:
                return path, key, f"CHECKSUM {etag} != {hex_md5}"
        else:
            # store's ETag isn't an MD5; size matched, so pass with a note
            LOGGER.debug("Non-MD5 ETag for %s, size-only verified", key)
        return None

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(check, paths):
            if result is not None:
                path, key, reason = result
                with lock:
                    problems["n"] += 1
                out.write(f"{reason}\t{key}\t{path}\n")

    LOGGER.info("Verified %d files, %d problems", len(paths), problems["n"])
    return problems["n"]


def _discover_prefixes(s3, bucket, root, depth, workers):
    def children(pfx):
        paginator = s3.get_paginator("list_objects_v2")
        out = []
        for page in paginator.paginate(Bucket=bucket, Prefix=pfx, Delimiter="/"):
            out += [c["Prefix"] for c in page.get("CommonPrefixes", [])]
        return out

    prefixes = [root]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in range(depth):
            expanded = []
            for p, kids in zip(prefixes, pool.map(children, prefixes)):
                expanded += kids if kids else [p]
            if expanded == prefixes:
                break
            prefixes = expanded
    return prefixes


def list_svgs(env, out, endpoint=ENDPOINT, bucket=BUCKET, workers=32, depth=2):
    """Write every SVG's bare URS in the bucket (one per line) to `out`.

    Parallelised by discovering {env}/URS/XX.. prefixes and paging them in a pool.
    Completeness is independent of depth: scan lists each prefix without a
    delimiter, so it recurses fully; depth only tunes fan-out for the skew.
    """
    s3 = client(endpoint, pool_connections=workers + 8)
    root = f"{env}/URS/"
    prefixes = _discover_prefixes(s3, bucket, root, depth, workers)
    LOGGER.info("Discovered %d prefixes under %s", len(prefixes), root)

    lock = threading.Lock()
    matched = {"n": 0}

    def scan(prefix):
        paginator = s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            urs = [
                name[: -len(SUFFIX)]
                for obj in page.get("Contents", [])
                for name in (obj["Key"].rsplit("/", 1)[-1],)
                if name.endswith(SUFFIX)
            ]
            if urs:
                with lock:
                    out.write("\n".join(urs) + "\n")
                    matched["n"] += len(urs)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for fut in as_completed([pool.submit(scan, p) for p in prefixes]):
            fut.result()
    LOGGER.info("Listed %d SVGs", matched["n"])
    return matched["n"]
