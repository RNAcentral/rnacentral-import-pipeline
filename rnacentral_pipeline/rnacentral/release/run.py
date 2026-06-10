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

LOGGER = logging.getLogger(__name__)

_BASE_CONNECT = {
    "keepalives": 1,
    "keepalives_idle": 60,
    "keepalives_interval": 10,
    "keepalives_count": 5,
}

# Conservative default: allows spilling to disk rather than OOM-killing the backend.
_CONNECT_DEFAULT = {**_BASE_CONNECT, "options": "-c statement_timeout=0 -c work_mem=64MB"}
# Higher memory only for DDL-heavy steps (index builds, partition exchange).
_CONNECT_HIGH_MEM = {**_BASE_CONNECT, "options": "-c statement_timeout=0 -c work_mem=256MB"}


def _connect(db_url, high_mem=False):
    return psycopg2.connect(db_url, **(_CONNECT_HIGH_MEM if high_mem else _CONNECT_DEFAULT))


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

TO_RELEASE = """
SELECT dbid, id
FROM rnacen.rnc_release
WHERE status = 'L'
ORDER BY id
"""

COUNT_QUERY = """
SELECT
    db.descr,
    count(distinct xref.upi)
from xref
join rnc_database db
on
    db.id = xref.dbid
where
    xref.deleted = 'N'
group by db.descr
"""

LOAD_COUNT_QUERY = """
SELECT
    load.database,
    count(distinct load.md5)
from load_rnacentral load
group by database
"""

PATCH_XREF_PARTITION_EXCHANGE_SQL = """
CREATE OR REPLACE FUNCTION rnc_load_xref.do_pel_exchange(p_in_db_id bigint)
RETURNS void
LANGUAGE plpgsql
SECURITY DEFINER
AS $function$
DECLARE
    v_parent_partition text;
    v_deleted_partition text := 'xref_p' || p_in_db_id || '_deleted';
    v_not_deleted_partition text := 'xref_p' || p_in_db_id || '_not_deleted';
    v_constraint_name text;
BEGIN
    --------------------------------------------------------
    --  drop any existing old partition tables
    --------------------------------------------------------
    execute 'drop table if exists xref_p' || p_in_db_id || '_deleted_old';
    execute 'drop table if exists xref_p' || p_in_db_id || '_not_deleted_old';

    --------------------------------------------------------
    --  rename indexes on current partition tables as old
    --------------------------------------------------------
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$id rename to xref_p' || p_in_db_id || '_deleted$id_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$dbid rename to xref_p' || p_in_db_id || '_deleted$dbid_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$ac rename to xref_p' || p_in_db_id || '_deleted$ac_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$created rename to xref_p' || p_in_db_id || '_deleted$created_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$last rename to xref_p' || p_in_db_id || '_deleted$last_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$upi rename to xref_p' || p_in_db_id || '_deleted$upi_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$taxid rename to xref_p' || p_in_db_id || '_deleted$taxid_old';

    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$id rename to xref_p' || p_in_db_id || '_not_deleted$id_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$dbid rename to xref_p' || p_in_db_id || '_not_deleted$dbid_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$ac rename to xref_p' || p_in_db_id || '_not_deleted$ac_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$created rename to xref_p' || p_in_db_id || '_not_deleted$created_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$last rename to xref_p' || p_in_db_id || '_not_deleted$last_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$upi rename to xref_p' || p_in_db_id || '_not_deleted$upi_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$taxid rename to xref_p' || p_in_db_id || '_not_deleted$taxid_old';

    --------------------------------------------------------
    --  rename FK constraints on current partition tables as old
    --------------------------------------------------------
    FOREACH v_constraint_name IN ARRAY ARRAY[
        'xref_p' || p_in_db_id || '_deleted_fk1',
        'xref_p' || p_in_db_id || '_deleted_fk2',
        'xref_p' || p_in_db_id || '_deleted_fk3',
        'xref_p' || p_in_db_id || '_deleted_fk4'
    ]
    LOOP
        IF EXISTS (
            SELECT 1
            FROM pg_constraint c
            JOIN pg_class r ON r.oid = c.conrelid
            JOIN pg_namespace n ON n.oid = r.relnamespace
            WHERE n.nspname = current_schema()
              AND r.relname = v_deleted_partition
              AND c.conname = v_constraint_name
        ) THEN
            execute format(
                'alter table %I rename constraint %I to %I',
                v_deleted_partition,
                v_constraint_name,
                v_constraint_name || '_old'
            );
        END IF;
    END LOOP;

    FOREACH v_constraint_name IN ARRAY ARRAY[
        'xref_p' || p_in_db_id || '_not_deleted_fk1',
        'xref_p' || p_in_db_id || '_not_deleted_fk2',
        'xref_p' || p_in_db_id || '_not_deleted_fk3',
        'xref_p' || p_in_db_id || '_not_deleted_fk4'
    ]
    LOOP
        IF EXISTS (
            SELECT 1
            FROM pg_constraint c
            JOIN pg_class r ON r.oid = c.conrelid
            JOIN pg_namespace n ON n.oid = r.relnamespace
            WHERE n.nspname = current_schema()
              AND r.relname = v_not_deleted_partition
              AND c.conname = v_constraint_name
        ) THEN
            execute format(
                'alter table %I rename constraint %I to %I',
                v_not_deleted_partition,
                v_constraint_name,
                v_constraint_name || '_old'
            );
        END IF;
    END LOOP;

    --------------------------------------------------------
    -- find the actual partition parent for the current dbid branch
    --------------------------------------------------------
    select inhparent::regclass::text
      into v_parent_partition
      from pg_inherits
     where inhrelid = to_regclass(v_deleted_partition);

    if v_parent_partition is null then
        raise exception 'Could not determine xref partition parent for dbid %', p_in_db_id;
    end if;

    --------------------------------------------------------
    -- detach current partitions before renaming them away
    --------------------------------------------------------
    execute format('alter table %s detach partition %I', v_parent_partition, v_deleted_partition);
    execute format('alter table %s detach partition %I', v_parent_partition, v_not_deleted_partition);

    --------------------------------------------------------
    --  Rename current partition tables as old
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename to xref_p' || p_in_db_id || '_deleted_old';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename to xref_p' || p_in_db_id || '_not_deleted_old';

    --------------------------------------------------------
    --  Rename pel tables according to naming conventions for partition tables of xref
    --------------------------------------------------------
    execute 'alter table xref_pel_deleted rename to xref_p' || p_in_db_id || '_deleted';
    execute 'alter table xref_pel_not_deleted rename to xref_p' || p_in_db_id || '_not_deleted';

    --------------------------------------------------------
    --  Add check constraints to accept only data complying with partitions definition
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_check' ||
            ' check ((dbid = ' || p_in_db_id || ') AND (deleted = ''Y''))';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_check' ||
            ' check ((dbid = ' || p_in_db_id || ') AND (deleted = ''N''))';

    --------------------------------------------------------
    --  Creating indexes on deleted partition
    --------------------------------------------------------
    execute 'create unique index xref_p' || p_in_db_id || '_deleted$id on xref_p' || p_in_db_id || '_deleted (id)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$dbid on xref_p' || p_in_db_id || '_deleted (dbid)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$ac on xref_p' || p_in_db_id || '_deleted (ac)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$created on xref_p' || p_in_db_id || '_deleted (created)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$last on xref_p' || p_in_db_id || '_deleted (last)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$upi on xref_p' || p_in_db_id || '_deleted (upi)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$taxid on xref_p' || p_in_db_id || '_deleted (taxid)';

    --------------------------------------------------------
    --  Creating FK constraints on deleted partition
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk1' ||
            ' foreign key(created) references rnc_release (id)';
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk2' ||
            ' foreign key(dbid) references rnc_database (id)';
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk3' ||
            ' foreign key(last) references rnc_release (id)';
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk4' ||
            ' foreign key(upi) references rna (upi)';

    --------------------------------------------------------
    --  Creating indexes on NOT deleted partition
    --------------------------------------------------------
    execute 'create unique index xref_p' || p_in_db_id || '_not_deleted$id on xref_p' || p_in_db_id || '_not_deleted (id)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$dbid on xref_p' || p_in_db_id || '_not_deleted (dbid)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$ac on xref_p' || p_in_db_id || '_not_deleted (ac)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$created on xref_p' || p_in_db_id || '_not_deleted (created)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$last on xref_p' || p_in_db_id || '_not_deleted (last)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$upi on xref_p' || p_in_db_id || '_not_deleted (upi)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$taxid on xref_p' || p_in_db_id || '_not_deleted (taxid)';

    --------------------------------------------------------
    --  Creating FK constraints on NOT deleted partition
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk1' ||
            ' foreign key(created) references rnc_release (id)';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk2' ||
            ' foreign key(dbid) references rnc_database (id)';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk3' ||
            ' foreign key(last) references rnc_release (id)';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk4' ||
            ' foreign key(upi) references rna (upi)';

    --------------------------------------------------------
    --  attach fresh partitions back to the correct parent
    --------------------------------------------------------
    execute format('alter table %s attach partition %I for values in (''Y'')', v_parent_partition, v_deleted_partition);
    execute format('alter table %s attach partition %I for values in (''N'')', v_parent_partition, v_not_deleted_partition);

    --------------------------------------------------------
    --   gathering statistics
    --------------------------------------------------------
    execute 'analyze xref_p' || p_in_db_id || '_deleted';
    execute 'analyze xref_p' || p_in_db_id || '_not_deleted';
END;
$function$;
"""


def run(db_url):
    """
    Run the release logic. Each step uses its own connection so a server-side
    crash on one long-running function doesn't abort the rest.
    """
    _run(db_url, PATCH_XREF_PARTITION_EXCHANGE_SQL, label="patch_xref_partition_exchange", high_mem=True)
    _run(db_url, "SELECT rnc_update.update_rnc_accessions()", label="update_rnc_accessions")
    _run(db_url, "SELECT rnc_update.update_literature_references()", label="update_literature_references")
    _run(db_url, CREATE_INDEX_SQL, label="create_index", high_mem=True)
    _run(db_url, "SELECT rnc_update.prepare_releases('F')", label="prepare_releases")

    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(TO_RELEASE)
            releases = cur.fetchall()

    for (dbid, rid) in releases:
        LOGGER.info("Executing release %i from database %i", rid, dbid)
        _run(db_url, "SELECT rnc_update.new_update_release(%s, %s)", params=(dbid, rid),
             label=f"new_update_release(dbid={dbid}, rid={rid})")


def check(limit_file, db_url, default_allowed_change=0.30):
    """
    Check the load tables for reasonable looking sequence counts.
    """

    limits = json.load(limit_file)
    cur_counts = {}
    new_counts = {}
    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(COUNT_QUERY)
            for (descr, raw_count) in cur.fetchall():
                cur_counts[descr] = float(raw_count)

            cur.execute(LOAD_COUNT_QUERY)
            for (descr, raw_count) in cur.fetchall():
                new_counts[descr] = float(raw_count)

    problems = False
    for name, previous in cur_counts.items():
        current = new_counts.get(name, default_allowed_change)
        change = (current - previous) / float(current)
        if change > limits.get(name, default_allowed_change):
            LOGGER.error("Database %s increased by %f", name, change)
            problems = True

    if problems:
        raise ValueError("Overly large changes with release")
