LOAD CSV
FROM ALL FILENAMES MATCHING ~<traveler-attempted.*csv$>
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
