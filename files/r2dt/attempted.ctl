LOAD CSV
FROM ALL FILENAMES MATCHING ~<r2dt-attempted.*csv$>
HAVING FIELDS (
  urs
)
INTO {{PGDATABASE}}?load_traveler_attempted
TARGET COLUMNS (
  urs
)

WITH
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_traveler_attempted;
$$,
$$
CREATE TABLE load_traveler_attempted (
  urs text primary key,
  r2dt_version text,
);
$$

AFTER LOAD DO
$$
INSERT INTO pipeline_tracking_traveler (
  urs,
  last_run,
  r2dt_version
) (
SELECT
  load.urs,
  NOW(),
  load.r2dt_version
FROM load_traveler_attempted load
) ON CONFLICT (urs) DO UPDATE
SET
  last_run = EXCLUDED.last_run
  r2dt_version = EXCLUDED.r2dt_version
;
$$
;
