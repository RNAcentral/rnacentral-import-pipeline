-- Copyright [2009-2014] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

set define off

create or replace PACKAGE RNC_LOGGING AS

  /* Package for logging RNAcentral operations. */

  PROCEDURE log_release_start(p_dbid NUMBER,
                              p_this_release NUMBER);

  PROCEDURE log_release_end(p_dbid NUMBER,
                            p_this_release NUMBER,
                            p_prev_release NUMBER := NULL);

END RNC_LOGGING;
/

create or replace PACKAGE BODY RNC_LOGGING AS


  /*
    Populate table RELEASE_STATS with basic stats about the staging table.
  */
  PROCEDURE log_release_start (
    p_dbid number,
    p_this_release number
  )
  IS
    PRAGMA AUTONOMOUS_TRANSACTION;
    l_dbid         NUMBER := p_dbid;

    l_this_release NUMBER:= p_this_release;
    l_sql_stmt     VARCHAR2 (4000);
  BEGIN

    l_sql_stmt := '
    MERGE INTO RELEASE_STATS s
    USING (
      SELECT
        :l_dbid dbid,
        :l_this_release this_release,
        SYSDATE start_time,
        COUNT (taxid) ff_taxid_count,
        COUNT (*) ff_loaded_rows -- all rows from the staging table

      FROM RNACEN.load_rnacentral
    ) q
    ON (s.this_release = q.this_release)
    -- reloading the same release
    WHEN MATCHED THEN UPDATE
      SET
      s.dbid = q.dbid,
      s.start_time = q.start_time,
      s.ff_taxid_nulls = (q.ff_loaded_rows - q.ff_taxid_count),
      s.ff_loaded_rows = q.ff_loaded_rows
    -- loading a new releaes
    WHEN NOT MATCHED THEN INSERT
      (

        dbid,
        this_release,
        start_time,
        ff_taxid_nulls,
        ff_loaded_rows
      )
      VALUES
      (
        q.dbid,
        q.this_release,
        q.start_time,
        (q.ff_loaded_rows - q.ff_taxid_count),
        q.ff_loaded_rows

      )'
    ;

    EXECUTE IMMEDIATE l_sql_stmt USING l_dbid, l_this_release;

    COMMIT;

  END log_release_start;

  /*

  */
  PROCEDURE log_release_end (

    p_dbid number,
    p_this_release NUMBER,
    p_prev_release number := NULL
  )
  IS
    PRAGMA AUTONOMOUS_TRANSACTION;
    l_dbid number := p_dbid;
    l_this_release NUMBER := p_this_release;
    l_prev_release number := NVL (NVL(p_prev_release, RNACEN.release.get_previous_release(l_dbid, l_this_release)), 0);
    l_sql_stmt varchar2 (4000);
  BEGIN

    MERGE INTO RELEASE_STATS s

    USING (
    SELECT
      l_dbid dbid,
      l_this_release this_release,
      l_prev_release prev_release,
      SYSDATE end_time,
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
        sum (retired_prev_releases) retired_prev_releases,

        sum (retired_this_release) retired_this_release,
        sum (retired_next_releases) retired_next_releases,
        sum (retired_total) retired_total,
        sum (case when version_i = 1 then created_w_predecessors else 0 end) created_w_predecessors_v_1,
        sum (case when version_i > 1 then created_w_predecessors else 0 end) created_w_predecessors_v_gt1,
        sum (created_w_predecessors) created_w_predecessors,
        sum (case when version_i = 1 then created_wo_predecessors else 0 end) created_wo_predecessors_v_1,
        sum (case when version_i > 1 then created_wo_predecessors else 0 end) created_wo_predecessors_v_gt1,
        sum (created_wo_predecessors) created_wo_predecessors,
        sum (active_created_prev_releases) active_created_prev_releases,
        sum (active_created_this_release) active_created_this_release,
        sum (active_created_next_releases) active_created_next_releases,
        sum (created_this_release) created_this_release,

        sum (active_updated_this_release) active_updated_this_release,
        sum (active_untouched_this_release) active_untouched_this_release,
        sum (active_total) active_total
          FROM (
            SELECT
              version_i,
              sum (
                CASE WHEN deleted = 'Y' and last < l_prev_release
                THEN 1
                ELSE 0
                END) retired_prev_releases,
              sum (
                CASE

                WHEN deleted = 'Y' AND LAST = l_prev_release
                THEN 1
                ELSE 0
                END) retired_this_release,
              sum (
                CASE
                WHEN deleted = 'Y' and last > l_prev_release
                THEN 1
                ELSE 0
                END) retired_next_releases,
              sum (
                CASE
                WHEN deleted = 'Y'

                then 1
                else 0
                END) retired_total,
              sum (
                CASE
                WHEN created = l_this_release AND EXISTS (
                  SELECT
                    1
                  FROM
                    xref p
                  WHERE
                    p.ac        = x.ac
                  and p.dbid    = x.dbid

                  and p.created < l_this_release)
                then 1
                else 0
                end) created_w_predecessors,
              sum (
                CASE
                WHEN created = l_this_release and not exists (
                  SELECT
                    1
                  FROM
                    xref p
                  WHERE
                    p.ac        = x.ac

                  and p.dbid    = x.dbid
                  and p.created < l_this_release)
                then 1
                else 0
                END) created_wo_predecessors,
              sum (
                CASE
                when deleted = 'N' and created < l_this_release
                then 1
                else 0
                END) active_created_prev_releases,
              sum (
                CASE

                when deleted = 'N' and created = l_this_release
                then 1
                else 0
                END) active_created_this_release,
              sum (
                case when deleted = 'N' and created > l_this_release
                then 1
                else 0
                end) active_created_next_releases,
              sum (
                case
                when created = l_this_release
                then 1

                else 0
                end) created_this_release,
              sum (
                case
                when deleted = 'N' and created != l_this_release and last = l_this_release
                then 1
                else 0
                end) active_updated_this_release,
              sum (
                case
                when deleted = 'N' and created != l_this_release and last != l_this_release
                then 1
                else 0

                end) active_untouched_this_release,
              sum (
                case
                WHEN deleted = 'N'
                then 1
                else 0
                END) active_total
        FROM xref x
        WHERE x.dbid = l_dbid
        GROUP BY version_i))
      ) q
    ON (s.this_release = q.this_release)
    WHEN MATCHED THEN UPDATE

    SET
      s.dbid = q.dbid,
      s.prev_release = q.prev_release,
      s.end_time = q.end_time,
      s.retired_prev_releases = q.retired_prev_releases,
      s.retired_this_release = q.retired_this_release,
      s.retired_next_releases = q.retired_next_releases,
      s.retired_total = q.retired_total,
      s.created_w_predecessors_v_1 = q.created_w_predecessors_v_1,
      s.created_w_predecessors_v_gt1 = q.created_w_predecessors_v_gt1,
      s.created_w_predecessors = q.created_w_predecessors,
      s.created_wo_predecessors_v_1 = q.created_wo_predecessors_v_1,
      s.created_wo_predecessors_v_gt1 = q.created_wo_predecessors_v_gt1,

      s.created_wo_predecessors = q.created_wo_predecessors,
      s.active_created_prev_releases = q.active_created_prev_releases,
      s.active_created_this_release = q.active_created_this_release,
      s.active_created_next_releases = q.active_created_next_releases,
      s.created_this_release = q.created_this_release,
      s.active_updated_this_release = q.active_updated_this_release,
      s.active_untouched_this_release = q.active_untouched_this_release,
      s.active_total = q.active_total
    WHEN NOT MATCHED THEN INSERT
    (
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
    VALUES
    (
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
    );

    COMMIT;

  END log_release_end;

END RNC_LOGGING;
/
set define on
