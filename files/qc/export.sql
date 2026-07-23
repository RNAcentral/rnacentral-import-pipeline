-- Export release summary: before->after active sequences per database and
-- rna_type, from rnc_rna_precomputed (before = not (re)precomputed this run, so
-- change is increase-only). CHECK if a row changed >5%. Vars: run_start.
SET statement_timeout = '600s';
\pset border 0
\pset footer off

\echo
\echo '### Sequences per database (before -> after this release; CHECK if >5% changed this run)'
SELECT database, before_sequences, after_sequences,
       round(100.0 * (after_sequences - before_sequences) / nullif(before_sequences, 0), 1) AS change_pct,
       CASE
         WHEN before_sequences = 0 AND after_sequences > 0 THEN 'CHECK: new (was 0)'
         WHEN (after_sequences - before_sequences) > 0.05 * before_sequences THEN 'CHECK: >5% change'
         ELSE 'ok'
       END AS status
FROM (
  SELECT db                                                                                AS database,
         count(*) FILTER (WHERE update_date < NULLIF(:'run_start', '')::timestamptz::date) AS before_sequences,
         count(*)                                                                          AS after_sequences
  FROM rnc_rna_precomputed p,
       unnest(string_to_array(p.databases, ',')) AS db
  WHERE p.is_active AND p.taxid IS NOT NULL
  GROUP BY db
) t
ORDER BY after_sequences DESC;

\echo
\echo '### Active sequences per rna_type (before -> after this release; CHECK if >5% changed this run)'
SELECT rna_type, before_sequences, after_sequences,
       round(100.0 * (after_sequences - before_sequences) / nullif(before_sequences, 0), 1) AS change_pct,
       CASE
         WHEN before_sequences = 0 AND after_sequences > 0 THEN 'CHECK: new (was 0)'
         WHEN (after_sequences - before_sequences) > 0.05 * before_sequences THEN 'CHECK: >5% change'
         ELSE 'ok'
       END AS status
FROM (
  SELECT coalesce(rna_type, '(none)')                                                      AS rna_type,
         count(*) FILTER (WHERE update_date < NULLIF(:'run_start', '')::timestamptz::date) AS before_sequences,
         count(*)                                                                          AS after_sequences
  FROM rnc_rna_precomputed
  WHERE is_active AND taxid IS NOT NULL
  GROUP BY 1
) t
ORDER BY after_sequences DESC;
