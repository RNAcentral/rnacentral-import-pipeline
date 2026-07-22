-- Genes QC: active / updated / new gene counts from rnc_genes (touched this run
-- = last_update = :release; flag if 0). Vars: release.
\pset pager off
\timing off
SET statement_timeout = '120s';

CREATE TEMP TABLE g AS
SELECT
  count(*) FILTER (WHERE is_active)                          AS active_genes,
  count(DISTINCT taxid) FILTER (WHERE is_active)             AS active_taxa,
  count(*) FILTER (WHERE last_update = NULLIF(:'release', '')::int)   AS updated_this_run,
  count(*) FILTER (WHERE first_release = NULLIF(:'release', '')::int) AS new_this_run
FROM rnc_genes;

\echo
\echo '### Genes — this run'
SELECT metric, value, status
FROM (
  SELECT 1 AS ord, 'active_genes' AS metric, active_genes::text AS value, 'ok' AS status FROM g
  UNION ALL
  SELECT 2, 'active_taxa', active_taxa::text, 'ok' FROM g
  UNION ALL
  SELECT 3, 'genes_updated_this_run', updated_this_run::text,
         CASE WHEN NULLIF(:'release', '') IS NOT NULL AND updated_this_run = 0 THEN 'CHECK' ELSE 'ok' END
  FROM g
  UNION ALL
  SELECT 4, 'new_genes_this_run', new_this_run::text, 'ok' FROM g
) m
ORDER BY ord;

\echo
\echo 'status  ok | CHECK: genes_updated_this_run=0 (genes step did not run / failed)'
