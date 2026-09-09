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

import json
import logging

import psycopg2

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
# Higher memory only for DDL-heavy steps (index builds, partition exchange).
_CONNECT_HIGH_MEM = {
    **_BASE_CONNECT,
    "options": "-c statement_timeout=0 -c work_mem=256MB",
}


def _connect(db_url, high_mem=False):
    return db.connect(db_url, **(_CONNECT_HIGH_MEM if high_mem else _CONNECT_DEFAULT))


def _run(db_url, sql, params=None, label="query", high_mem=False):
    with _connect(db_url, high_mem=high_mem) as conn:
        conn.autocommit = True
        with conn.cursor() as cur:
            LOGGER.info("Running %s", label)
            cur.execute(sql, params)


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

# Restricted to databases actually staged in this run. An abandoned load leaves a
# pending release behind, and without this filter that stale release gets loaded
# instead of the staged one.
TO_RELEASE = """
SELECT r.dbid, r.id
FROM rnacen.rnc_release r
WHERE r.status = 'L'
AND r.dbid IN (
    SELECT DISTINCT d.id
    FROM rnacen.load_rnacentral_all l
    JOIN rnacen.rnc_database d ON l.database = d.descr
)
ORDER BY r.id
"""

LOAD_COUNT_QUERY = """
SELECT
    load.database,
    count(distinct load.md5)
from load_rnacentral load
group by database
"""

# Databases staged for this run, so COUNT_QUERY below can filter to their
# partitions instead of aggregating every xref_pN partition regardless of
# what this run actually loaded.
LOADED_DBIDS_QUERY = """
SELECT DISTINCT d.id
FROM load_rnacentral l
JOIN rnc_database d ON d.descr = l.database
"""

# xref.dbid = ANY(%s) prunes to just this run's partitions - previously an
# unfiltered GROUP BY over all ~43 databases' worth of xref, every release,
# even though only the loaded ones are ever compared against LOAD_COUNT_QUERY.
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
    and xref.dbid = ANY(%s)
group by db.descr
"""


def run(db_url):
    """
    Run the release logic. Each step uses its own connection so a server-side
    crash on one long-running function doesn't abort the rest.
    """
    # Deploy any changed database functions from database_functions/ before the
    # release logic runs. Replaces the per-function CREATE OR REPLACE patches that
    # used to live here inline.
    functions.apply(db_url)
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
    _run(db_url, "SELECT rnc_update.prepare_releases('F')", label="prepare_releases")

    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(TO_RELEASE)
            releases = cur.fetchall()

    for (dbid, rid) in releases:
        LOGGER.info("Executing release %i from database %i", rid, dbid)
        _run(
            db_url,
            "SELECT rnc_update.new_update_release(%s, %s)",
            params=(dbid, rid),
            label=f"new_update_release(dbid={dbid}, rid={rid})",
        )

    # do_pel_exchange adds each partition's upi->rna foreign key (fk4) NOT VALID to keep the
    # full-partition validation scan off the load's critical path. Validate them now, after
    # the per-database loads have committed: VALIDATE CONSTRAINT takes only
    # ShareUpdateExclusiveLock (does not block readers/writers), marks the constraint valid,
    # and surfaces any (by-construction vanishingly unlikely) violation. The data is already
    # committed at this point, so this is detection, not a pre-commit gate.
    for (dbid, rid) in releases:
        for suffix in ("deleted", "not_deleted"):
            _run(
                db_url,
                f"ALTER TABLE xref_p{dbid}_{suffix} "
                f"VALIDATE CONSTRAINT xref_p{dbid}_{suffix}_fk4",
                label=f"validate fk4 xref_p{dbid}_{suffix}",
            )

    # Verify xref primary key uniqueness for each database loaded this run.
    # do_checks scopes its scan to the dbid's own partitions against the rest
    # of xref, rather than aggregating the whole table, so a per-database call
    # is cheap.
    for (dbid, rid) in releases:
        _run(
            db_url,
            "SELECT rnc_load_xref.do_checks(%s::bigint)",
            params=(dbid,),
            label=f"do_checks(dbid={dbid})",
            high_mem=True,
        )


def check(limit_file, db_url, default_allowed_change=0.30):
    """
    Check the load tables for reasonable looking sequence counts.
    """

    limits = json.load(limit_file)
    cur_counts = {}
    new_counts = {}
    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(LOAD_COUNT_QUERY)
            for (descr, raw_count) in cur.fetchall():
                new_counts[descr] = float(raw_count)

            cur.execute(LOADED_DBIDS_QUERY)
            dbids = [dbid for (dbid,) in cur.fetchall()]

            cur.execute(COUNT_QUERY, (dbids,))
            for (descr, raw_count) in cur.fetchall():
                cur_counts[descr] = float(raw_count)

    problems = False
    for name, previous in cur_counts.items():
        current = new_counts.get(name, default_allowed_change)
        change = (current - previous) / float(current)
        if change > limits.get(name, default_allowed_change):
            LOGGER.error("Database %s increased by %f", name, change)
            problems = True

    if problems:
        raise ValueError("Overly large changes with release")
