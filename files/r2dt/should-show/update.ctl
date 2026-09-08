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
  inferred_should_show = load.should_show
FROM load_secondary_should_show load
WHERE
  layout.urs = load.urs
;
$$,
$$
DROP TABLE load_secondary_should_show;
$$
;
