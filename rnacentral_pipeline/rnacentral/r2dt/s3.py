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
import csv
import gzip
import hashlib
import logging
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3
import pyarrow.parquet as pq
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

# Sized for moving an object's bytes.
READ_TIMEOUT = 120
# Listing is different: a page is a few tens of KB, so anything slower than this
# means the store is struggling, not that the response is large. Waiting the full
# READ_TIMEOUT there just parks a worker for two minutes before the retry that
# was always going to be needed -- with every worker doing it, throughput
# collapses while the connection pressure that caused it stays pinned.
LIST_READ_TIMEOUT = 30


def client(endpoint=ENDPOINT, pool_connections=64, read_timeout=READ_TIMEOUT):
    """A path-style boto3 S3 client for the Weka endpoint.

    Honours S3_KEY/S3_SECRET (as update-svg.sh did) first, else falls back to
    boto3's normal chain (~/.aws/credentials, AWS_PROFILE, instance role, ...).

    `read_timeout` defaults to a value sized for transferring an object. Listing
    passes something much shorter: see LIST_READ_TIMEOUT.
    """
    kwargs = dict(
        endpoint_url=endpoint,
        config=Config(
            s3={"addressing_style": "path"},
            max_pool_connections=pool_connections,
            connect_timeout=30,
            read_timeout=read_timeout,
            # One attempt, because _with_retries already gives each object five.
            # Nesting the two budgets multiplied them: a wedged key cost
            # 5 x 5 x 120s ~= 50 minutes before the task gave up. Adaptive mode
            # still throttles client-side when the store pushes back.
            retries={"max_attempts": 1, "mode": "adaptive"},
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


def _already_stored(s3, bucket, key, hex_md5, size):
    """True when `key` already holds exactly these bytes.

    A read timeout on a PUT usually means the object landed and only the
    response was lost, so re-PUTting it times out the same way -- which is how
    one object failed all MAX_ATTEMPTS and took a whole publish_layout chunk
    down with it (1/3124 uploads failed). A HEAD is cheap; if the bytes are
    there the upload is done. Any error here means "can't tell", so the caller
    falls back to retrying.
    """
    try:
        head = s3.head_object(Bucket=bucket, Key=key)
    except Exception:  # noqa: BLE001 - unreachable/404/denied all mean "retry"
        return False
    if head.get("ContentLength") != size:
        return False
    etag = head.get("ETag", "").strip('"')
    if _looks_like_md5(etag):
        return etag == hex_md5
    # Store's ETag isn't an MD5 (same caveat as verify()); size matched.
    return True


def _with_retries(fn, key, attempts=MAX_ATTEMPTS, sleep=None, settled=None):
    """Call `fn`, retrying transient errors with exponential backoff.

    `settled`, if given, is consulted after a failure and should report whether
    the work actually landed despite the error; when it does we stop and return
    None rather than retrying.
    """
    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except Exception as err:  # noqa: BLE001 - classified by _is_transient
            if settled is not None and _is_transient(err) and settled():
                LOGGER.warning(
                    "%s reported %s but landed anyway; not retrying", key, err
                )
                return None
            if attempt == attempts or not _is_transient(err):
                raise
            LOGGER.warning(
                "Retrying %s after %s (attempt %d/%d)", key, err, attempt, attempts
            )
            # Resolved per call, not as a default argument: a default would bind
            # time.sleep at import and quietly ignore a monkeypatched sleep.
            (sleep or time.sleep)(BACKOFF * 2 ** (attempt - 1))


def upload(
    file_list,
    env,
    endpoint=ENDPOINT,
    bucket=BUCKET,
    workers=32,
    allow_failures=0,
    failure_list=None,
):
    """Upload every gzipped SVG named in `file_list` to S3, in parallel.

    Each PUT sends Content-MD5 (server-side integrity) and we assert the returned
    ETag matches (client-side). Raises RuntimeError if ANY object failed, so the
    caller/pipeline sees the failure rather than a silent partial upload.

    `allow_failures` relaxes that for stray objects the store will not accept --
    a single wedged key (every PUT and HEAD to it timing out) otherwise costs the
    whole batch and, after the retries, the run. Failures up to the cap are
    logged and written to `failure_list`, one URS per line, so the gap is
    recorded and can be reconciled with verify-s3 rather than silently lost.
    Anything above the cap still raises: an outage must not slip through as
    tolerated breakage.
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
            settled=lambda: _already_stored(s3, bucket, key, hex_md5, len(data)),
        )
        if resp is None:  # the PUT landed, only its response was lost
            return key
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
        if failure_list is not None:
            Path(failure_list).write_text(
                "".join(f"{urs_of(path)}\t{err}\n" for path, err in failures)
            )
        summary = (
            f"{len(failures)}/{len(paths)} uploads failed; first: "
            f"{failures[0][0]}: {failures[0][1]}"
        )
        if len(failures) > allow_failures:
            raise RuntimeError(summary)
        LOGGER.error("%s -- tolerated (allow_failures=%d)", summary, allow_failures)
    return done["n"]


def _looks_like_md5(etag):
    return len(etag) == 32 and all(c in "0123456789abcdef" for c in etag.lower())


def failed_urs(failure_list):
    """The URS ids in a --failure-list file (one per line, URS first).

    `#` lines are the run headers the published, appended-to copy carries, so
    that file can be read back here as well as by a person.
    """
    path = Path(failure_list)
    if not path.is_file():
        return set()
    return {
        line.split("\t", 1)[0].strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def _drop_from_csv(path, failed):
    rows = list(csv.reader(path.open()))
    keep = [row for row in rows if not row or row[0] not in failed]
    if len(keep) == len(rows):
        return 0
    with path.open("w", newline="") as out:
        csv.writer(out).writerows(keep)
    return len(rows) - len(keep)


def _drop_from_parquet(path, failed):
    table = pq.read_table(path)
    if "urs" not in table.column_names:
        return 0
    kept = table.filter([u not in failed for u in table.column("urs").to_pylist()])
    if kept.num_rows == table.num_rows:
        return 0
    pq.write_table(kept, path)
    return table.num_rows - kept.num_rows


def drop_failed(failure_list, paths, max_failures=None):
    """
    Strip rows for URS whose SVG never reached S3 from the files about to load.

    Both sides must lose the row: left in attempted*.parquet it reaches
    pipeline_tracking_traveler, which files/r2dt/setup.sql excludes from future
    runs, so the sequence would never be drawn again. Above max_failures this
    raises -- that many missing objects is an outage, not a stray wedged key.
    """
    failed = failed_urs(failure_list)
    if not failed:
        return 0

    if max_failures is not None and len(failed) > max_failures:
        raise RuntimeError(
            f"{len(failed)} uploads failed across the run, over the "
            f"{max_failures} allowed: {', '.join(sorted(failed))}"
        )

    dropped = 0
    for name in paths:
        path = Path(name)
        if path.suffix == ".parquet":
            dropped += _drop_from_parquet(path, failed)
        else:
            dropped += _drop_from_csv(path, failed)

    LOGGER.warning(
        "Dropped %d row(s) for %d URS whose SVG did not upload: %s",
        dropped,
        len(failed),
        ", ".join(sorted(failed)),
    )
    return dropped


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


def _pages(s3, bucket, prefix, delimiter=None):
    """Yield every listing page under `prefix`, retrying a failed page in place.

    Deliberately not paginator.paginate: listing was the one S3 path here with no
    retry at all, while client() sets max_attempts=1 on the assumption that
    _with_retries covers everything -- so a single read timeout on one page threw
    away a whole multi-million-object inventory run.

    A retry has to resume from the last continuation token rather than restart
    the prefix: the caller streams keys straight out, and sync() COPYs that into
    a `urs text PRIMARY KEY` table, which a re-emitted key would fail.
    """
    token = None
    pages = 0
    while True:
        kwargs = {"Bucket": bucket, "Prefix": prefix}
        if delimiter:
            kwargs["Delimiter"] = delimiter
        if token:
            kwargs["ContinuationToken"] = token
        page = _with_retries(lambda: s3.list_objects_v2(**kwargs), prefix)
        pages += 1
        LOGGER.debug(
            "%s page %d: %d keys, truncated=%s, token=%s",
            prefix,
            pages,
            len(page.get("Contents", [])),
            page.get("IsTruncated"),
            bool(page.get("NextContinuationToken")),
        )
        yield page

        truncated = page.get("IsTruncated")
        token = page.get("NextContinuationToken")
        if not truncated:
            return
        # Truncated but no token to continue from: the store is telling us there
        # is more and refusing to say where. Stopping here would silently drop
        # the rest of the prefix and still report a clean total -- exactly the
        # class of silent divergence this module exists to prevent. Fail instead.
        if not token:
            raise RuntimeError(
                f"{prefix} page {pages} is truncated but carries no "
                "NextContinuationToken; listing would silently be short"
            )


def _discover_prefixes(s3, bucket, root, depth, workers):
    def children(pfx):
        out = []
        for page in _pages(s3, bucket, pfx, delimiter="/"):
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
    s3 = client(endpoint, pool_connections=workers + 8, read_timeout=LIST_READ_TIMEOUT)
    root = f"{env}/URS/"
    prefixes = _discover_prefixes(s3, bucket, root, depth, workers)
    LOGGER.info("Discovered %d prefixes under %s", len(prefixes), root)

    lock = threading.Lock()
    matched = {"n": 0}
    done = {"n": 0}

    def scan(prefix):
        for page in _pages(s3, bucket, prefix):
            urs = [
                name[: -len(SUFFIX)]
                for obj in page.get("Contents", [])
                for name in (obj["Key"].rsplit("/", 1)[-1],)
                if name.endswith(SUFFIX)
            ]
            if urs:
                with lock:
                    out.write("\n".join(urs) + "\n")
                    # Flushed per page: a listing run takes tens of minutes and
                    # the default 8KB buffer holds back the only evidence it is
                    # progressing, which reads as a hang. One flush per ~1000
                    # keys is free next to the request that fetched them.
                    out.flush()
                    matched["n"] += len(urs)
                    done["n"] += 1
                    if done["n"] % 100 == 0:
                        LOGGER.info("%d pages, %d SVGs listed", done["n"], matched["n"])

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for fut in as_completed([pool.submit(scan, p) for p in prefixes]):
            fut.result()
    LOGGER.info("Listed %d SVGs", matched["n"])
    return matched["n"]


# The should-show set, matching the partial index r2dt_results_urs_show_idx and
# every other consumer (precompute, search-export): a manual assigned_should_show
# overrides the model's inferred_should_show.
SHOULD_SHOW_SQL = """
SELECT urs
FROM r2dt_results
WHERE coalesce(assigned_should_show, inferred_should_show)
GROUP BY urs
"""

DELETE_BATCH = 1000  # S3 DeleteObjects caps at 1000 keys per call


def _copy_out(conn, query, path):
    """Stream a single-column query straight to `path` via COPY, counting rows.

    COPY rather than a cursor loop because these sets run to tens of millions of
    rows: psycopg2 writes the stream out without ever building a Python list.
    """
    with conn.cursor() as cur, open(path, "w") as out:
        cur.copy_expert(f"COPY ({query}) TO STDOUT", out)
        return cur.rowcount


def sync(
    env,
    db_url,
    missing_list,
    orphan_list,
    delete=False,
    endpoint=ENDPOINT,
    bucket=BUCKET,
    workers=32,
    depth=2,
):
    """Reconcile the bucket against the DB should-show set.

    Writes two files: `missing_list` (should show but no SVG in S3 -- the list to
    feed back through r2dt) and `orphan_list` (in S3 but no longer should show).
    Only removes the orphans from S3 when `delete` is true; the default is a
    dry run, because a mis-set --env would otherwise wipe the wrong environment.

    The set difference happens in Postgres, not in Python: at ~36M objects a pair
    of in-memory sets would cost several GB, whereas an indexed anti-join is flat.
    """
    from rnacentral_pipeline import db

    inventory = Path(missing_list).with_suffix(".s3-inventory.txt")
    with open(inventory, "w") as out:
        in_s3 = list_svgs(
            env, out, endpoint=endpoint, bucket=bucket, workers=workers, depth=depth
        )

    with db.connection(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute("DROP TABLE IF EXISTS s3_svg_urs")
            cur.execute("CREATE UNLOGGED TABLE s3_svg_urs (urs text PRIMARY KEY)")
            with open(inventory) as handle:
                cur.copy_expert("COPY s3_svg_urs (urs) FROM STDIN", handle)
            cur.execute("ANALYZE s3_svg_urs")

        missing = _copy_out(
            conn,
            f"SELECT d.urs FROM ({SHOULD_SHOW_SQL}) d "
            "LEFT JOIN s3_svg_urs s USING (urs) WHERE s.urs IS NULL",
            missing_list,
        )
        orphan = _copy_out(
            conn,
            f"SELECT s.urs FROM s3_svg_urs s "
            f"LEFT JOIN ({SHOULD_SHOW_SQL}) d USING (urs) WHERE d.urs IS NULL",
            orphan_list,
        )

        with conn.cursor() as cur:
            cur.execute("DROP TABLE s3_svg_urs")

    LOGGER.info(
        "in_s3=%d missing=%d orphan=%d (%s)",
        in_s3,
        missing,
        orphan,
        "deleting orphans" if delete else "dry run, nothing deleted",
    )
    if delete and orphan:
        delete_orphans(
            orphan_list, env, endpoint=endpoint, bucket=bucket, workers=workers
        )
    return {"in_s3": in_s3, "missing": missing, "orphan": orphan}


def _batched_lines(path, size):
    """Yield lists of up to `size` stripped, non-empty lines. Streams the file."""
    batch = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            batch.append(line)
            if len(batch) == size:
                yield batch
                batch = []
    if batch:
        yield batch


def delete_orphans(
    orphan_list, env, endpoint=ENDPOINT, bucket=BUCKET, workers=32, batch=DELETE_BATCH
):
    """Delete every URS named in `orphan_list` from the bucket.

    Batched into DeleteObjects calls and retried per batch. Raises if the store
    reports any per-key error, so a partial prune is never mistaken for a clean
    one. Only ever called behind sync(delete=True) / --delete.
    """
    s3 = client(endpoint, pool_connections=workers + 8)
    lock = threading.Lock()
    done = {"n": 0}
    failures = []

    def delete_batch(urs_batch):
        keys = [{"Key": s3_key(u, env)} for u in urs_batch]
        resp = _with_retries(
            lambda: s3.delete_objects(Bucket=bucket, Delete={"Objects": keys}),
            keys[0]["Key"],
        )
        errors = (resp or {}).get("Errors") or []
        if errors:
            raise RuntimeError(
                f"{len(errors)} keys failed to delete; first: {errors[0]}"
            )
        return len(keys)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(delete_batch, b) for b in _batched_lines(orphan_list, batch)
        ]
        for fut in as_completed(futures):
            try:
                count = fut.result()
            except Exception as err:  # noqa: BLE001 - collect, report all at end
                with lock:
                    failures.append(err)
                LOGGER.error("Failed to delete a batch: %s", err)
            else:
                with lock:
                    done["n"] += count

    LOGGER.info("Deleted %d objects from %s", done["n"], bucket)
    if failures:
        raise RuntimeError(
            f"{len(failures)} delete batches failed; first: {failures[0]}"
        )
    return done["n"]


def local_path(directory, urs):
    """Where a URS lands on disk: the S3 key minus the env prefix, ungzipped.

    Sharded the same way as the bucket because a flat directory of tens of
    millions of files is unusable for most filesystems and for tar.
    """
    key = s3_key(urs, "_")  # env prefix stripped below; only the shards matter
    relative = key.split("/", 1)[1][: -len(SUFFIX)] + ".svg"
    return Path(directory) / relative


def download(
    directory,
    env,
    urs_list=None,
    endpoint=ENDPOINT,
    bucket=BUCKET,
    workers=32,
    depth=2,
    manifest="manifest.tsv",
):
    """Download every SVG to `directory`, ungzipped, for bulk distribution.

    Pulls the URS from `urs_list` if given (e.g. the reconciled should-show list),
    otherwise from a fresh bucket listing. The list is streamed, never held in
    memory. Files are written decompressed so the eventual `tar czf` actually
    compresses -- re-gzipping .svg.gz gains nothing.

    Resumable: an existing non-empty destination is skipped, so re-running after
    an interruption is cheap. That also means a resumed run appends to the
    manifest, which may then hold duplicate rows for re-fetched URS.

    Raises RuntimeError if any object failed, matching upload()'s contract.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    if urs_list is None:
        urs_list = directory / "s3-inventory.txt"
        with open(urs_list, "w") as out:
            list_svgs(
                env, out, endpoint=endpoint, bucket=bucket, workers=workers, depth=depth
            )

    s3 = client(endpoint, pool_connections=workers + 8)
    lock = threading.Lock()
    counts = {"done": 0, "skipped": 0}
    failures = []
    manifest_path = directory / manifest

    def fetch(urs):
        dest = local_path(directory, urs)
        if dest.exists() and dest.stat().st_size:
            with lock:
                counts["skipped"] += 1
            return None
        key = s3_key(urs, env)
        resp = _with_retries(lambda: s3.get_object(Bucket=bucket, Key=key), key)
        svg = gzip.decompress(resp["Body"].read())
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(svg)
        hex_md5, _ = _md5(svg)
        return f"{urs}\t{dest.relative_to(directory)}\t{len(svg)}\t{hex_md5}\n"

    with open(manifest_path, "a") as manifest_fh, ThreadPoolExecutor(
        max_workers=workers
    ) as pool:
        for batch in _batched_lines(urs_list, 10000):
            futures = {pool.submit(fetch, u): u for u in batch}
            for fut in as_completed(futures):
                urs = futures[fut]
                try:
                    row = fut.result()
                except Exception as err:  # noqa: BLE001 - collect, report at end
                    with lock:
                        failures.append((urs, err))
                    LOGGER.error("Failed to download %s: %s", urs, err)
                else:
                    if row is None:
                        continue
                    with lock:
                        manifest_fh.write(row)
                        counts["done"] += 1

    LOGGER.info(
        "Downloaded %d, skipped %d already present, %d failed",
        counts["done"],
        counts["skipped"],
        len(failures),
    )
    if failures:
        raise RuntimeError(
            f"{len(failures)} downloads failed; first: "
            f"{failures[0][0]}: {failures[0][1]}"
        )
    return counts["done"]
