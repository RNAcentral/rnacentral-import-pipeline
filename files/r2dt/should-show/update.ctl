LOAD CSV
FROM ALL FILENAMES MATCHING ~<r2dt-should-show.*.csv$>
HAVING FIELDS (
  urs,
  should_show
) INTO {{PGDATABASE}}?load_secondary_should_show
TARGET COLUMNS (
  urs,
  should_show
)

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_secondary_should_show;
$$,
$$
CREATE UNLOGGED TABLE load_secondary_should_show (
  urs text NOT NULL,
  should_show bool NOT NULL
);
$$

WITH
  FIELDS ESCAPED BY double-quote,
  FIELDS TERMINATED BY ','

AFTER LOAD DO
$$
UPDATE r2dt_results layout
SET
  -- The classifier must not re-show a force-fit that load.ctl stored hidden.
  inferred_should_show = load.should_show
    AND coalesce(rna.len::numeric / nullif(m.model_length, 0) < 1.25, true)
FROM load_secondary_should_show load, rna, r2dt_models m
WHERE
  layout.urs = load.urs
  AND rna.urs = layout.urs
  AND m.id = layout.model_id
;
$$,
$$
DROP TABLE load_secondary_should_show;
$$
;
