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

with pred as (
    -- Per-ac predecessor flag for this dbid, computed in a single pass.
    -- Replaces two per-row correlated EXISTS(...) subqueries against xref that
    -- each probed xref(ac, dbid, created) once per outer row, turning an
    -- O(rows * lookups) nested loop into one hash aggregate + hash join.
    -- An ac "has a predecessor" iff any of its xref rows (for this dbid) was
    -- created before this release, i.e. the earliest created < l_this_release.
    select ac, (min(created) < l_this_release) as has_predecessor
    from xref
    where dbid = l_dbid
    group by ac
),
new_stats as (
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
                WHEN created = l_this_release AND hp.has_predecessor
                then 1
                else 0
                end) created_w_predecessors,
              sum(
                CASE
                WHEN created = l_this_release and not hp.has_predecessor
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
        -- pred is built from the same xref/dbid rows, so every x.ac is present;
        -- the join always matches (has_predecessor is never NULL for these rows).
        JOIN pred hp ON hp.ac = x.ac
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
$function$

