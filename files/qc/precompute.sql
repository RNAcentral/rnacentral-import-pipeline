-- Precompute QC: active-sequence counts + flags, from rnc_rna_precomputed and
-- qa_status. active = is_active AND taxid IS NOT NULL; precomputed this run =
-- update_date >= run start. Vars: run_start, qa_issue_pct_max (default 50).
SET statement_timeout = '240s';
\pset border 0
\pset footer off

\if :{?qa_issue_pct_max}
\else
\set qa_issue_pct_max 50
\endif

CREATE TEMP TABLE pc AS
SELECT
  count(*) FILTER (WHERE taxid IS NOT NULL) AS active_sequences,
  count(*) FILTER (WHERE taxid IS NOT NULL
                     AND update_date >= NULLIF(:'run_start', '')::timestamptz::date) AS precomputed_this_run,
  count(*) FILTER (WHERE taxid IS NOT NULL AND rna_type IS NULL) AS active_missing_rna_type
FROM rnc_rna_precomputed
WHERE is_active;

CREATE TEMP TABLE qa AS
SELECT count(*) AS total, count(*) FILTER (WHERE has_issue) AS with_issue
FROM qa_status;

\echo
SELECT metric, value, status
FROM (
  SELECT 1 AS ord, 'active_sequences' AS metric, active_sequences::text AS value, 'ok' AS status
  FROM pc
  UNION ALL
  SELECT 2, 'precomputed_this_run', precomputed_this_run::text,
         CASE WHEN NULLIF(:'run_start', '') IS NOT NULL AND precomputed_this_run = 0 THEN 'CHECK: precompute did not run' ELSE 'ok' END
  FROM pc
  UNION ALL
  SELECT 3, 'active_missing_rna_type', active_missing_rna_type::text,
         CASE WHEN active_missing_rna_type > 0 THEN 'CHECK: rna_type assignment incomplete' ELSE 'ok' END
  FROM pc
  UNION ALL
  SELECT 4, 'qa_sequences_with_issues', round(100.0 * with_issue / nullif(total, 0), 1)::text || '%',
         CASE WHEN 100.0 * with_issue / nullif(total, 0) > :qa_issue_pct_max THEN 'CHECK: above sanity bound' ELSE 'ok' END
  FROM qa
) m
ORDER BY ord;
