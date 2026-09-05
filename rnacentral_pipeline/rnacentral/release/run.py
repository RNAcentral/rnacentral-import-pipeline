# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import datetime
import json
import logging
import select
import time

import psycopg2
import psycopg2.extensions
from retry import retry

from rnacentral_pipeline import db
from rnacentral_pipeline.rnacentral.release import functions

LOGGER = logging.getLogger(__name__)

_BASE_CONNECT = {
    "keepalives": 1,
    "keepalives_idle": 60,
    "keepalives_interval": 10,
    "keepalives_count": 5,
}

# Conservative default: allows spilling to disk rather than OOM-killing the backend.
_CONNECT_DEFAULT = {
    **_BASE_CONNECT,
    "options": "-c statement_timeout=0 -c work_mem=64MB",
}
# Higher memory for the steps that sort or hash at table scale. maintenance_work_mem
# must be set explicitly: CREATE INDEX and VALIDATE CONSTRAINT ignore work_mem, so
# setting work_mem alone left those steps on the server default. Kept under the 2GB
# populate_precompute.sql already asks of this server, where index builds have OOM'd.
_CONNECT_HIGH_MEM = {
    **_BASE_CONNECT,
    "options": "-c statement_timeout=0 -c work_mem=1GB -c maintenance_work_mem=1GB",
}


# Seconds to block in select() between polls, and how often to log a liveness
# line while a statement is still running. The release functions can go hours
# between notices, so without the heartbeat a running step is indistinguishable
# from a hung one.
_POLL_TIMEOUT = 5.0
_HEARTBEAT_SECONDS = 15 * 60


def _connect(db_url, high_mem=False):
    return db.connect(db_url, **(_CONNECT_HIGH_MEM if high_mem else _CONNECT_DEFAULT))


def _elapsed(started):
    return datetime.timedelta(seconds=round(time.monotonic() - started))


def _log_notices(conn, label):
    """
    Drain server-side messages (the RAISE NOTICE progress lines in the release
    functions) into the log. Returns True if anything was logged.
    """
    logged = False
    notices = getattr(conn, "notices", None)
    while notices:
        notice = notices.pop(0).strip()
        if notice:
            LOGGER.info("%s | %s", label, notice)
            logged = True
    return logged


def _wait(conn, label):
    """
    Drive an async connection to completion, logging notices as they arrive.

    psycopg2 only converts pending notices into conn.notices when the connection
    is polled, and a synchronous execute() polls nothing until the statement has
    finished -- so a multi-hour function's progress notices all arrive at once,
    at the end, which is useless. Polling by hand surfaces each one as the server
    sends it.
    """
    started = time.monotonic()
    last_beat = started
    while True:
        state = conn.poll()
        if _log_notices(conn, label):
            # A notice is itself proof of life; don't follow it with a heartbeat
            # saying the same thing.
            last_beat = time.monotonic()
        if state == psycopg2.extensions.POLL_OK:
            return
        if state == psycopg2.extensions.POLL_READ:
            select.select([conn.fileno()], [], [], _POLL_TIMEOUT)
        elif state == psycopg2.extensions.POLL_WRITE:
            select.select([], [conn.fileno()], [], _POLL_TIMEOUT)
        else:
            raise psycopg2.OperationalError("Unexpected poll state %s" % state)
        now = time.monotonic()
        if now - last_beat >= _HEARTBEAT_SECONDS:
            last_beat = now
            LOGGER.info("STILL  %s (%s)", label, _elapsed(started))


@retry(
    psycopg2.OperationalError, tries=6, delay=5, backoff=2, max_delay=120, jitter=(0, 5)
)
def _connect_async(db_url, high_mem=False, label="connect"):
    """
    Open an async connection. Async mode has no implicit transaction, so each
    statement commits on its own -- the same behaviour the old autocommit
    connection had. The handshake finishes inside _wait, not inside connect(),
    which is why the retry lives here as well as in db.connect.
    """
    options = _CONNECT_HIGH_MEM if high_mem else _CONNECT_DEFAULT
    conn = db.connect(db_url, async_=1, **options)
    try:
        _wait(conn, label)
    except Exception:
        conn.close()
        raise
    return conn


def _run(db_url, sql, params=None, label="query", high_mem=False):
    LOGGER.info("START  %s", label)
    started = time.monotonic()
    conn = _connect_async(db_url, high_mem=high_mem, label=label)
    try:
        with conn.cursor() as cur:
            cur.execute(sql, params)
            _wait(conn, label)
    finally:
        _log_notices(conn, label)
        conn.close()
    LOGGER.info("DONE   %s (%s)", label, _elapsed(started))


CREATE_INDEX_SQL = """
CREATE INDEX IF NOT EXISTS load_rnacentral_all$database
ON rnacen.load_rnacentral_all(database)
"""

# Index on load_md5_new_sequences(in_md5), the join key used by set_comparable_prot_upi
# and store_new_sequences. load_md5_new_sequences is only TRUNCATEd (not dropped) during
# a load, so a one-time index here survives every per-database iteration. (The functions
# now join rather than probe per-row, so the planner may hash-join instead; the index is
# kept as a cheap safety net and is harmless if unused.)
LOAD_MD5_INDEX_SQL = """
CREATE INDEX IF NOT EXISTS load_md5_new_sequences$in_md5
ON rnacen.load_md5_new_sequences(in_md5)
"""

TO_RELEASE = """
SELECT dbid, id
FROM rnacen.rnc_release
WHERE status = 'L'
ORDER BY id
"""

COUNT_QUERY = """
SELECT
    db.descr,
    count(distinct xref.urs)
from xref
join rnc_database db
on
    db.id = xref.dbid
where
    xref.deleted = 'N'
    and xref.dbid in ({dbids})
group by db.descr
"""

LOAD_COUNT_QUERY = """
SELECT
    load.database,
    count(distinct load.md5)
from load_rnacentral load
group by database
"""


def run(db_url, force_full=False):
    """
    Run the release logic. Each step uses its own connection so a server-side
    crash on one long-running function doesn't abort the rest.

    By default each database's release type is chosen automatically: the first
    load of a database is FULL (it bootstraps the xref partition), every later
    load is INCREMENTAL (only new/changed rows are touched). Pass force_full=True
    to force a FULL release for every database -- e.g. after a schema change.
    """
    run_started = time.monotonic()
    # Deploy any changed database functions from database_functions/ before the
    # release logic runs. Replaces the per-function CREATE OR REPLACE patches that
    # used to live here inline.
    LOGGER.info("START  apply_functions")
    started = time.monotonic()
    functions.apply(db_url)
    LOGGER.info("DONE   apply_functions (%s)", _elapsed(started))
    _run(
        db_url,
        "SELECT rnc_update.update_rnc_accessions()",
        label="update_rnc_accessions",
    )
    _run(
        db_url,
        "SELECT rnc_update.update_literature_references()",
        label="update_literature_references",
    )
    _run(db_url, CREATE_INDEX_SQL, label="create_index", high_mem=True)
    _run(db_url, LOAD_MD5_INDEX_SQL, label="create_load_md5_index", high_mem=True)
    # 'A' = auto (per-database FULL/INCREMENTAL from history); 'F' forces FULL.
    _run(
        db_url,
        "SELECT rnc_update.prepare_releases(%s)",
        params=("F" if force_full else "A",),
        label="prepare_releases",
    )

    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(TO_RELEASE)
            releases = cur.fetchall()

    LOGGER.info(
        "%i database(s) to release: %s",
        len(releases),
        ", ".join(str(dbid) for (dbid, _rid) in releases),
    )
    for (n, (dbid, rid)) in enumerate(releases, start=1):
        LOGGER.info("Executing release %i from database %i", rid, dbid)
        _run(
            db_url,
            "SELECT rnc_update.new_update_release(%s, %s)",
            params=(dbid, rid),
            label=f"new_update_release(dbid={dbid}, rid={rid}) [{n}/{len(releases)}]",
        )

    # do_pel_exchange adds each partition's urs->rna foreign key (fk4) NOT VALID to keep the
    # full-partition validation scan off the load's critical path. Validate them now, after
    # the per-database loads have committed: VALIDATE CONSTRAINT takes only
    # ShareUpdateExclusiveLock (does not block readers/writers), marks the constraint valid,
    # and surfaces any (by-construction vanishingly unlikely) violation. The data is already
    # committed at this point, so this is detection, not a pre-commit gate.
    for (dbid, rid) in releases:
        for suffix in ("deleted", "not_deleted"):
            # Validation joins the whole partition against rna, so on ENA this is
            # the most spill-prone statement in the release tail.
            _run(
                db_url,
                f"ALTER TABLE xref_p{dbid}_{suffix} "
                f"VALIDATE CONSTRAINT xref_p{dbid}_{suffix}_fk4",
                label=f"validate fk4 xref_p{dbid}_{suffix}",
                high_mem=True,
            )

    # Verify xref primary key uniqueness once, after all databases are loaded,
    # rather than once per database inside load_xref. The check is global (it
    # ignores its argument), so a single run covers every partition.
    if releases:
        _run(
            db_url,
            "SELECT rnc_load_xref.do_checks(NULL::bigint)",
            label="do_checks (once, post-loop)",
            high_mem=True,
        )

    LOGGER.info("DONE   release run (%s total)", _elapsed(run_started))


def check(limit_file, db_url, default_allowed_change=0.30):
    """
    Check the load tables for reasonable looking sequence counts.

    Two-sided: the same per-database limit (from limit_file, else
    default_allowed_change) applies to both growth and shrink. The shrink side
    matters because incremental/full releases retire any active xref ABSENT from
    the load, so a parse that silently comes up far short would delete real rows by
    absence. Smaller swings are expected and left for the QC report to surface.

    Delta-parsed databases are exempt from the shrink check: their load table holds
    only new/changed rows (deletions come from an explicit list, not absence), so a
    loaded count far below the active count is expected, not a red flag.
    """

    limits = json.load(limit_file)
    cur_counts = {}
    new_counts = {}
    delta_dbs = set()
    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(LOAD_COUNT_QUERY)
            for (descr, raw_count) in cur.fetchall():
                new_counts[descr] = float(raw_count)

            # Only loaded databases can change; scope the xref count to their
            # partitions (dbids spliced as literals so the partitions prune).
            if new_counts:
                cur.execute(
                    "SELECT id FROM rnc_database WHERE descr = ANY(%s)",
                    (list(new_counts),),
                )
                dbids = ",".join(str(int(row[0])) for row in cur)
                if dbids:
                    cur.execute(COUNT_QUERY.format(dbids=dbids))
                    for (descr, raw_count) in cur.fetchall():
                        cur_counts[descr] = float(raw_count)

            # Databases loaded via delta parsing (an import manifest exists for them);
            # exempt from the shrink check below. Guard the lookup: the manifest table
            # may not exist on older schemas.
            cur.execute("SELECT to_regclass('rnacen.pipeline_tracking_import')")
            if cur.fetchone()[0] is not None:
                cur.execute("SELECT database FROM rnacen.pipeline_tracking_import")
                delta_dbs = {row[0] for row in cur}

    problems = False
    for name, previous in cur_counts.items():
        current = new_counts.get(name, 0.0)
        if not current:
            continue
        limit = limits.get(name, default_allowed_change)
        increase = (current - previous) / float(current)
        if increase > limit:
            LOGGER.error("Database %s increased by %f", name, increase)
            problems = True
        elif name not in delta_dbs and previous:
            decrease = (previous - current) / float(previous)
            if decrease > limit:
                LOGGER.error(
                    "Database %s decreased by %f (loaded %d distinct sequences vs %d active)",
                    name,
                    decrease,
                    int(current),
                    int(previous),
                )
                problems = True

    if problems:
        raise ValueError("Overly large changes with release")
