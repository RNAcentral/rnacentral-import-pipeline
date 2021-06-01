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
  urs text primary key
);
$$

AFTER LOAD DO
$$
INSERT INTO pipeline_tracking_traveler (
  urs,
  last_run
) (
SELECT
  load.urs,
  NOW()
FROM load_traveler_attempted load
) ON CONFLICT (urs) DO UPDATE
SET 
  last_run = EXCLUDED.last_run
;
$$
;
