-- Import QC: per-database rows imported / created / retired this release, from
-- release_stats (correct for ENA). Flags NOT IMPORTED (import_date < run start)
-- and "data disappeared" (>retired_pct retired AND carried_over dropped).
-- Vars: run_dbs, run_start, retired_pct (default 0.05).
\pset pager off
\timing off
SET statement_timeout = '60s';

\if :{?retired_pct}
\else
\set retired_pct 0.05
\endif

WITH hist AS (
  SELECT dbid, this_release, prev_release, end_time,
         ff_loaded_rows, created_this_release, retired_this_release,
         active_created_prev_releases, active_total,
         lag(active_created_prev_releases) OVER w AS prev_carried_over,
         lag(active_total)                 OVER w AS prev_active_total
  FROM release_stats
  WINDOW w AS (PARTITION BY dbid ORDER BY this_release)
),
latest AS (
  SELECT DISTINCT ON (dbid) *
  FROM hist
  ORDER BY dbid, this_release DESC
)
SELECT d.display_name              AS database,
       l.prev_release,
       l.this_release              AS current_release,
       l.end_time::date            AS import_date,
       l.ff_loaded_rows            AS rows_imported,
       l.created_this_release      AS created,
       l.retired_this_release      AS retired,
       CASE
         WHEN l.end_time < NULLIF(:'run_start', '')::timestamptz
           THEN 'NOT IMPORTED'
         WHEN l.prev_active_total > 0
              AND l.retired_this_release > :retired_pct * l.prev_active_total
              AND l.active_created_prev_releases < l.prev_carried_over
           THEN 'CHECK: data disappeared'
         ELSE 'ok'
       END                         AS status
FROM rnc_database d
JOIN latest l ON l.dbid = d.id
WHERE lower(replace(replace(d.descr, '_', ''), ' ', '')) = ANY (string_to_array(:'run_dbs', ','))
   OR lower(replace(replace(d.display_name, '_', ''), ' ', '')) = ANY (string_to_array(:'run_dbs', ','))
ORDER BY d.display_name;

\echo
\echo 'status  ok = imported cleanly | NOT IMPORTED = no new release this run (import_date < run start) | CHECK: data disappeared = >5% retired AND carried_over fell vs last release (possible partial/failed import)'
