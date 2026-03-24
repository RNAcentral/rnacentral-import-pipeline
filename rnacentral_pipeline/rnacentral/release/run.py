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
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk1 to xref_p' || p_in_db_id || '_deleted_fk1_old';
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk2 to xref_p' || p_in_db_id || '_deleted_fk2_old';
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk3 to xref_p' || p_in_db_id || '_deleted_fk3_old';
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk4 to xref_p' || p_in_db_id || '_deleted_fk4_old';

    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk1 to xref_p' || p_in_db_id || '_not_deleted_fk1_old';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk2 to xref_p' || p_in_db_id || '_not_deleted_fk2_old';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk3 to xref_p' || p_in_db_id || '_not_deleted_fk3_old';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk4 to xref_p' || p_in_db_id || '_not_deleted_fk4_old';

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


def patch_xref_partition_exchange(cursor):
    cursor.execute(PATCH_XREF_PARTITION_EXCHANGE_SQL)


def run(db_url):
    """
    Run the release logic. Basically this will run the select commands that are
    needed to update all the tables and then run the release update logic.
    """

    with psycopg2.connect(db_url) as conn:
        cursor = conn.cursor()
        cursor.execute("SET work_mem TO '256MB'")
        patch_xref_partition_exchange(cursor)
        cursor.execute("SELECT rnc_update.update_rnc_accessions()")
        cursor.execute("SELECT rnc_update.update_literature_references()")
        cursor.execute(CREATE_INDEX_SQL)
        cursor.execute("SELECT rnc_update.prepare_releases('F')")
        cursor.execute(TO_RELEASE)
        for (dbid, rid) in cursor.fetchall():
            LOGGER.info("Executing release %i from database %i", rid, dbid)
            cursor.execute("SELECT rnc_update.new_update_release(%s, %s)", (dbid, rid))
            conn.commit()


def check(limit_file, db_url, default_allowed_change=0.30):
    """
    Check the load tables for reasonable looking sequence counts.
    """

    limits = json.load(limit_file)
    cur_counts = {}
    new_counts = {}
    with psycopg2.connect(db_url) as conn:
        cursor = conn.cursor()
        cursor.execute(COUNT_QUERY)
        for (descr, raw_count) in cursor.fetchall():
            cur_counts[descr] = float(raw_count)

        cursor.execute(LOAD_COUNT_QUERY)
        for (descr, raw_count) in cursor.fetchall():
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
