LOAD CSV
FROM ALL FILENAMES MATCHING ~<data.*.csv$>
HAVING FIELDS (
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
) INTO {{PGDATABASE}}?load_secondary
TARGET COLUMNS (
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
)

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_secondary;
$$,
$$
CREATE UNLOGGED TABLE load_secondary (
  urs text primary key,
  model_id int,
  secondary_structure text,
  overlap_count int,
  basepair_count int,
  model_start int,
  model_stop int,
  sequence_start int,
  sequence_stop int,
  sequence_coverage float,
  inferred_should_show bool
);
$$

WITH
    batch rows = 300,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

AFTER LOAD DO
$$
-- The website reads should_show directly, so a force-fit (sequence far longer
-- than its model) has to be stored hidden, not just filtered out on export.
UPDATE load_secondary load
SET inferred_should_show = false
FROM rna, r2dt_models m
WHERE rna.urs = load.urs
  AND m.id = load.model_id
  AND rna.len::numeric / nullif(m.model_length, 0) >= 1.25;
$$,
$$
INSERT INTO r2dt_results (
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
) (
SELECT
    urs,
    model_id,
    secondary_structure,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    inferred_should_show
FROM load_secondary
) ON CONFLICT (urs) DO UPDATE
SET
    model_id = EXCLUDED.model_id,
    secondary_structure = EXCLUDED.secondary_structure,
    overlap_count = EXCLUDED.overlap_count,
    basepair_count = EXCLUDED.basepair_count,
    model_start = EXCLUDED.model_start,
    model_stop = EXCLUDED.model_stop,
    sequence_start = EXCLUDED.sequence_start,
    sequence_stop = EXCLUDED.sequence_stop,
    sequence_coverage = EXCLUDED.sequence_coverage,
    inferred_should_show = EXCLUDED.inferred_should_show
;
$$,
$$
DROP TABLE load_secondary
$$
;
