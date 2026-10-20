-- Analyze QC. -v snapshot=1 emits the before-snapshot CSV; otherwise reports
-- sequences run per analysis. Vars: snapshot, run_start.
SET statement_timeout = '180s';

\if :{?snapshot}

-- Before-snapshot: CSV of untracked result-table counts.
\pset tuples_only on
\pset format unaligned
\pset fieldsep ','
SELECT 'cpat',     count(*) FROM cpat_results
UNION ALL
SELECT 'stopfree', count(*) FROM stopfree_results
UNION ALL
SELECT 'tcode',    count(*) FROM tcode_results;

\else

\pset border 0
\pset footer off

-- After-report: tracked analyses via last_run, untracked via before/after diff.
CREATE TEMP TABLE before_counts (analysis text, before_count bigint);
\copy before_counts FROM 'analyze-before.csv' WITH (FORMAT csv)

\echo
WITH tracked AS (
  SELECT 'rfam scan' AS analysis, 1 AS ord,
         count(*) FILTER (WHERE last_run >= NULLIF(:'run_start', '')::timestamptz) AS sequences_run,
         count(*) AS running_total
  FROM pipeline_tracking_qa_scan
  UNION ALL
  SELECT 'r2dt', 2,
         count(*) FILTER (WHERE last_run >= NULLIF(:'run_start', '')::timestamptz), count(*)
  FROM pipeline_tracking_traveler
  UNION ALL
  SELECT 'genome mapping', 3,
         count(*) FILTER (WHERE last_run >= NULLIF(:'run_start', '')::timestamptz), count(*)
  FROM pipeline_tracking_genome_mapping
),
untracked AS (
  SELECT ac.analysis, ac.ord,
         ac.after_count - coalesce(bc.before_count, 0) AS sequences_run,
         ac.after_count AS running_total
  FROM (
    SELECT 'cpat'     AS analysis, 4 AS ord, count(*) AS after_count FROM cpat_results
    UNION ALL SELECT 'tcode',    5, count(*) FROM tcode_results
    UNION ALL SELECT 'stopfree', 6, count(*) FROM stopfree_results
  ) ac
  LEFT JOIN before_counts bc ON bc.analysis = ac.analysis
)
SELECT analysis, sequences_run, running_total,
       CASE
         WHEN NULLIF(:'run_start', '') IS NOT NULL AND sequences_run <= 0
           THEN 'CHECK: nothing processed (step skipped/failed?)'
         ELSE 'ok'
       END AS status
FROM (
  SELECT analysis, ord, sequences_run, running_total FROM tracked
  UNION ALL
  SELECT analysis, ord, sequences_run, running_total FROM untracked
) x
ORDER BY ord;

\endif
