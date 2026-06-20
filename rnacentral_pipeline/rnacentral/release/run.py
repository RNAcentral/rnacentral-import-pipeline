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

# rnc_load_rna.set_comparable_prot_upi runs a correlated subquery that probes
# load_md5_new_sequences by in_md5 once per load_retro_tmp row. The table has no
# index, so that is a per-row sequential scan. load_md5_new_sequences is only
# TRUNCATEd (not dropped) during a load, so a one-time index here survives every
# per-database iteration.
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


# do_checks runs `select id, count(*) from xref group by id having count(*) > 1`,
# a full scan + global aggregate of the *entire* xref table. The deployed
# load_xref calls it once per database, so a single release re-checks the whole
# table ~once per loaded database. It is logically a global invariant, so it only
# needs to run once. This patch removes the per-database call; run() invokes
# do_checks a single time after the per-database loop instead.
PATCH_LOAD_XREF_SQL = """
CREATE OR REPLACE FUNCTION rnc_load_xref.load_xref(p_previous_release bigint, p_in_dbid bigint)
RETURNS void
LANGUAGE plpgsql
SECURITY DEFINER
AS $function$
BEGIN
    perform rnc_load_xref.prepare_pel_tables();
    perform rnc_load_xref.populate_pel_tables1(p_previous_release);
    perform rnc_load_xref.populate_pel_tables2();
    perform rnc_load_xref.load_upi_max_versions_table(p_in_dbid);
    perform rnc_load_xref.load_max_versions_table();
    perform rnc_load_xref.populate_pel_tables3(p_in_dbid);
    perform rnc_load_xref.populate_pel_tables4(p_in_dbid, p_previous_release);
    perform rnc_load_xref.do_pel_exchange(p_in_dbid);
    -- do_checks intentionally NOT called here; run() calls it once after the
    -- per-database loop (it scans the whole xref table regardless of dbid).
END;
$function$;
"""


# The xref joins in populate_pel_tables1/2 and load_upi_max_versions_table filter
# the (LIST-by-dbid partitioned) xref table on dbid as a *join key*, never a
# constant, so the planner cannot prune and Appends across all ~59 partitions for
# every outer row. Supplying dbid as a value the planner can prune on collapses
# this to a single-partition scan (verified: "Subplans Removed: 58").
#
# All three additions are logically redundant: load_retro_tmp only ever holds
# rows for the one database being loaded (load_retro_tmp_table truncates it and
# inserts p_in_dbid as in_dbid for every row), so they cannot change which rows
# match -- they only enable partition pruning.
PATCH_POPULATE_PEL_TABLES1_SQL = """
CREATE OR REPLACE FUNCTION rnc_load_xref.populate_pel_tables1(v_previous_release bigint)
RETURNS void
LANGUAGE plpgsql
SECURITY DEFINER
AS $function$
DECLARE
    -- load_retro_tmp is single-dbid, so this is the dbid for the whole table.
    v_dbid bigint := (SELECT in_dbid FROM load_retro_tmp LIMIT 1);
BEGIN

with l_data AS
(
  select x.dbid, x.created, x.upi, x.version_i, x.timestamp,
         x.userstamp, x.ac, x.version, case
           when (x.last < l.in_load_release and x.upi = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null))) then l.in_load_release
           else coalesce(v_previous_release,x.last)
         end as last, case
           when (x.last < l.in_load_release and x.upi = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null))) then 'N'
           else 'Y'
         end as deleted,
         case
           when (x.last < l.in_load_release and x.upi = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null))) then coalesce(l.in_taxid,x.taxid)
           else x.taxid
         end taxid
  from load_retro_tmp l,
       xref x
  where x.ac = l.in_ac
  and   x.dbid = l.in_dbid
  and   x.dbid = v_dbid   -- partition pruning hint; redundant (load_retro_tmp is single-dbid)
  and   l.comparable_prot_upi is not null
  and   x.deleted = 'N'
),
ins1 AS
(
  insert into xref_pel_deleted
  (
    dbid,
    created,
    upi,
    version_i,
    timestamp,
    userstamp,
    ac,
    VERSION,
    last,
    deleted,
    taxid
  )
  select dbid, created, upi, version_i, timestamp,
         userstamp, ac, VERSION, last, deleted,
         taxid
  from l_data
  where deleted = 'Y'
)
insert into xref_pel_not_deleted
(
  dbid,
  created,
  upi,
  version_i,
  timestamp,
  userstamp,
  ac,
  version,
  last,
  deleted,
  taxid
)
select dbid, created, upi, version_i, timestamp,
       userstamp, ac, version, last, deleted,
       taxid
from l_data
where deleted = 'N';

END;
$function$;
"""

PATCH_POPULATE_PEL_TABLES2_SQL = """
CREATE OR REPLACE FUNCTION rnc_load_xref.populate_pel_tables2()
RETURNS void
LANGUAGE plpgsql
SECURITY DEFINER
AS $function$
DECLARE
    -- load_retro_tmp is single-dbid, so this is the dbid for the whole table.
    v_dbid bigint := (SELECT in_dbid FROM load_retro_tmp LIMIT 1);
BEGIN

  insert into xref_pel_not_deleted
  (
    dbid,
    created,
    upi,
    version_i,
    timestamp,
    userstamp,
    ac,
    version,
    last,
    deleted,
    taxid
  )
  select l.in_dbid, l.in_load_release created, l.comparable_prot_upi upi,
         case
           when (x.upi != l.comparable_prot_upi) then x.version_i + 1
           else x.version_i
         end,
         clock_timestamp() as "timestamp",
         user userstamp, l.in_ac ac, l.in_version as "version", l.in_load_release as "last", 'N' deleted,
         in_taxid taxid
  from load_retro_tmp l, xref x
  where x.ac = l.in_ac
  and   x.dbid = l.in_dbid
  and   x.dbid = v_dbid   -- partition pruning hint; redundant (load_retro_tmp is single-dbid)
  and   l.comparable_prot_upi is not null
  and   not -- this condition differentiates this procedure from populate_pel_tables1
        (x.last < l.in_load_release and x.upi = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null)))
  and   x.deleted = 'N';

END;
$function$;
"""

# Here dbid is already available as the function argument p_in_dbid, and the WHERE
# already enforces L.IN_DBID = p_in_dbid, so referencing p_in_dbid directly on the
# xref side is exact -- and gives the planner a value to prune on.
PATCH_LOAD_UPI_MAX_VERSIONS_SQL = """
CREATE OR REPLACE FUNCTION rnc_load_xref.load_upi_max_versions_table(p_in_dbid bigint)
RETURNS void
LANGUAGE plpgsql
SECURITY DEFINER
AS $function$
BEGIN

    EXECUTE 'TRUNCATE TABLE load_upi_max_versions';

    INSERT INTO load_upi_max_versions
      SELECT DISTINCT
        L.IN_AC,
        L.IN_DBID,
        MAX(coalesce(PREVIOUS_XREF.VERSION_I, 0)) MAX_VERSION_I, -- VERSION_I set to 0 at first
        PREVIOUS_XREF.UPI
      FROM rnacen.xref previous_xref RIGHT OUTER JOIN load_retro_tmp l ON (PREVIOUS_XREF.DBID = p_in_dbid AND PREVIOUS_XREF.AC = L.IN_AC)
     WHERE L.COMPARABLE_PROT_UPI IS NOT NULL AND L.IN_DBID = p_in_dbid   -- outer join, left side can be NULL
     -- outer join, left side can be NULL
       AND NOT EXISTS (
            SELECT 1
            FROM RNACEN.XREF X
            WHERE X.AC      = L.IN_AC
            AND X.DBID    = p_in_dbid
            AND X.DELETED = 'N'
          )
      GROUP BY
        L.IN_AC,
        L.IN_DBID,
        PREVIOUS_XREF.UPI;

    --COMMIT;

    execute 'analyze load_upi_max_versions';

  END;
$function$;
"""


# rnc_logging.log_release_end_atx (called per-database at the end of load_release via
# log_release_end) aggregates over xref WHERE x.dbid = l_dbid -- which prunes -- but its
# two correlated EXISTS subqueries filter the inner xref scan with p.dbid = x.dbid, a
# correlated column rather than a constant. That inner scan therefore Appends across all
# ~59 partitions for every outer row, observed live as a 46-minute, all-partition scan on
# a 630k-row database. Since the outer WHERE already forces x.dbid = l_dbid, replacing
# p.dbid = x.dbid with p.dbid = l_dbid is exact and lets the inner scan prune. Same family
# of fix as the populate_pel_tables/load_upi_max_versions patches above; this function lives
# in rnc_logging, which is not part of the database_functions dump we first analysed.
PATCH_LOG_RELEASE_END_ATX_SQL = """
CREATE OR REPLACE FUNCTION rnc_logging.log_release_end_atx(p_dbid bigint, p_this_release bigint, p_prev_release bigint DEFAULT NULL::bigint)
RETURNS void
LANGUAGE plpgsql
SECURITY DEFINER
AS $function$
DECLARE
    l_dbid RNACEN.rnc_release.dbid%TYPE := p_dbid;
    l_this_release RNACEN.rnc_release.id%TYPE := p_this_release;
    l_prev_release RNACEN.rnc_release.id%TYPE := coalesce(coalesce(p_prev_release, release.get_previous_release(l_dbid, l_this_release)), 0);
    l_sql_stmt varchar(4000);

BEGIN

with new_stats as (
    SELECT
      l_dbid dbid,
      l_this_release this_release,
      l_prev_release prev_release,
      clock_timestamp() end_time,
      retired_prev_releases,
      retired_this_release,
      retired_next_releases,
      retired_total,
      created_w_predecessors_v_1,
      created_w_predecessors_v_gt1,
      created_w_predecessors,
      created_wo_predecessors_v_1,
      created_wo_predecessors_v_gt1,
      created_wo_predecessors,
      active_created_prev_releases,
      active_created_this_release,
      active_created_next_releases,
      created_this_release,
      active_updated_this_release,
      active_untouched_this_release,
      active_total
    FROM (
      SELECT
        sum(retired_prev_releases) retired_prev_releases,
        sum(retired_this_release) retired_this_release,
        sum(retired_next_releases) retired_next_releases,
        sum(retired_total) retired_total,
        sum(case when version_i = 1 then created_w_predecessors else 0 end) created_w_predecessors_v_1,
        sum(case when version_i > 1 then created_w_predecessors else 0 end) created_w_predecessors_v_gt1,
        sum(created_w_predecessors) created_w_predecessors,
        sum(case when version_i = 1 then created_wo_predecessors else 0 end) created_wo_predecessors_v_1,
        sum(case when version_i > 1 then created_wo_predecessors else 0 end) created_wo_predecessors_v_gt1,
        sum(created_wo_predecessors) created_wo_predecessors,
        sum(active_created_prev_releases) active_created_prev_releases,
        sum(active_created_this_release) active_created_this_release,
        sum(active_created_next_releases) active_created_next_releases,
        sum(created_this_release) created_this_release,
        sum(active_updated_this_release) active_updated_this_release,
        sum(active_untouched_this_release) active_untouched_this_release,
        sum(active_total) active_total
          FROM (
            SELECT
              version_i,
              sum(
                CASE WHEN deleted = 'Y' and last < l_prev_release
                THEN 1
                ELSE 0
                END) retired_prev_releases,
              sum(
                CASE
                WHEN deleted = 'Y' AND LAST = l_prev_release
                THEN 1
                ELSE 0
                END) retired_this_release,
              sum(
                CASE
                WHEN deleted = 'Y' and last > l_prev_release
                THEN 1
                ELSE 0
                END) retired_next_releases,
              sum(
                CASE
                WHEN deleted = 'Y'
                then 1
                else 0
                END) retired_total,
              sum(
                CASE
                WHEN created = l_this_release AND EXISTS (
                  SELECT
                    1
                  FROM
                    xref p
                  WHERE
                    p.ac        = x.ac
                  and p.dbid    = l_dbid
                  and p.created < l_this_release)
                then 1
                else 0
                end) created_w_predecessors,
              sum(
                CASE
                WHEN created = l_this_release and not exists (
                  SELECT
                    1
                  FROM
                    xref p
                  WHERE
                    p.ac        = x.ac
                  and p.dbid    = l_dbid
                  and p.created < l_this_release)
                then 1
                else 0
                END) created_wo_predecessors,
              sum(
                CASE
                when deleted = 'N' and created < l_this_release
                then 1
                else 0
                END) active_created_prev_releases,
              sum(
                CASE
                when deleted = 'N' and created = l_this_release
                then 1
                else 0
                END) active_created_this_release,
              sum(
                case when deleted = 'N' and created > l_this_release
                then 1
                else 0
                end) active_created_next_releases,
              sum(
                case
                when created = l_this_release
                then 1
                else 0
                end) created_this_release,
              sum(
                case
                when deleted = 'N' and created != l_this_release and last = l_this_release
                then 1
                else 0
                end) active_updated_this_release,
              sum(
                case
                when deleted = 'N' and created != l_this_release and last != l_this_release
                then 1
                else 0
                end) active_untouched_this_release,
              sum(
                case
                WHEN deleted = 'N'
                then 1
                else 0
                END) active_total
        FROM xref as x
        WHERE x.dbid = l_dbid
        GROUP BY version_i) alias33
      ) q
    )
    insert INTO RELEASE_STATS as s (
      dbid,
      this_release,
      prev_release,
      end_time,
      retired_prev_releases,
      retired_this_release,
      retired_next_releases,
      retired_total,
      created_w_predecessors_v_1,
      created_w_predecessors_v_gt1,
      created_w_predecessors,
      created_wo_predecessors_v_1,
      created_wo_predecessors_v_gt1,
      created_wo_predecessors,
      active_created_prev_releases,
      active_created_this_release,
      active_created_next_releases,
      created_this_release,
      active_updated_this_release,
      active_untouched_this_release,
      active_total
    )
  select
       q.dbid,
       q.this_release,
       q.prev_release,
       q.end_time,
       q.retired_prev_releases,
       q.retired_this_release,
       q.retired_next_releases,
       q.retired_total,
       q.created_w_predecessors_v_1,
       q.created_w_predecessors_v_gt1,
       q.created_w_predecessors,
       q.created_wo_predecessors_v_1,
       q.created_wo_predecessors_v_gt1,
       q.created_wo_predecessors,
       q.active_created_prev_releases,
       q.active_created_this_release,
       q.active_created_next_releases,
       q.created_this_release,
       q.active_updated_this_release,
       q.active_untouched_this_release,
       q.active_total
  from new_stats as q
    on conflict (this_release) do update
    SET
      dbid = excluded.dbid,
      prev_release = excluded.prev_release,
      end_time = excluded.end_time,
      retired_prev_releases = excluded.retired_prev_releases,
      retired_this_release = excluded.retired_this_release,
      retired_next_releases = excluded.retired_next_releases,
      retired_total = excluded.retired_total,
      created_w_predecessors_v_1 = excluded.created_w_predecessors_v_1,
      created_w_predecessors_v_gt1 = excluded.created_w_predecessors_v_gt1,
      created_w_predecessors = excluded.created_w_predecessors,
      created_wo_predecessors_v_1 = excluded.created_wo_predecessors_v_1,
      created_wo_predecessors_v_gt1 = excluded.created_wo_predecessors_v_gt1,
      created_wo_predecessors = excluded.created_wo_predecessors,
      active_created_prev_releases = excluded.active_created_prev_releases,
      active_created_this_release = excluded.active_created_this_release,
      active_created_next_releases = excluded.active_created_next_releases,
      created_this_release = excluded.created_this_release,
      active_updated_this_release = excluded.active_updated_this_release,
      active_untouched_this_release = excluded.active_untouched_this_release,
      active_total = excluded.active_total;

  END;
$function$;
"""


def run(db_url):
    """
    Run the release logic. Each step uses its own connection so a server-side
    crash on one long-running function doesn't abort the rest.
    """
    _run(db_url, PATCH_XREF_PARTITION_EXCHANGE_SQL, label="patch_xref_partition_exchange", high_mem=True)
    _run(db_url, PATCH_LOAD_XREF_SQL, label="patch_load_xref", high_mem=True)
    _run(db_url, PATCH_POPULATE_PEL_TABLES1_SQL, label="patch_populate_pel_tables1", high_mem=True)
    _run(db_url, PATCH_POPULATE_PEL_TABLES2_SQL, label="patch_populate_pel_tables2", high_mem=True)
    _run(db_url, PATCH_LOAD_UPI_MAX_VERSIONS_SQL, label="patch_load_upi_max_versions", high_mem=True)
    _run(db_url, PATCH_LOG_RELEASE_END_ATX_SQL, label="patch_log_release_end_atx", high_mem=True)
    _run(db_url, "SELECT rnc_update.update_rnc_accessions()", label="update_rnc_accessions")
    _run(db_url, "SELECT rnc_update.update_literature_references()", label="update_literature_references")
    _run(db_url, CREATE_INDEX_SQL, label="create_index", high_mem=True)
    _run(db_url, LOAD_MD5_INDEX_SQL, label="create_load_md5_index", high_mem=True)
    _run(db_url, "SELECT rnc_update.prepare_releases('F')", label="prepare_releases")

    with _connect(db_url) as conn:
        with conn.cursor() as cur:
            cur.execute(TO_RELEASE)
            releases = cur.fetchall()

    for (dbid, rid) in releases:
        LOGGER.info("Executing release %i from database %i", rid, dbid)
        _run(db_url, "SELECT rnc_update.new_update_release(%s, %s)", params=(dbid, rid),
             label=f"new_update_release(dbid={dbid}, rid={rid})")

    # Verify xref primary key uniqueness once, after all databases are loaded,
    # rather than once per database inside load_xref. The check is global (it
    # ignores its argument), so a single run covers every partition.
    if releases:
        _run(db_url, "SELECT rnc_load_xref.do_checks(NULL::bigint)",
             label="do_checks (once, post-loop)", high_mem=True)


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
