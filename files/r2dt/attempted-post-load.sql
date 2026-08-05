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
  last_run = EXCLUDED.last_run,
  r2dt_version = EXCLUDED.r2dt_version;

-- Clear the staging table now that it is merged. Without this it is never
-- emptied between runs: it reached 1.6 GB (86k live rows behind a 1.5 GB
-- primary key) in production before this was added.
TRUNCATE load_traveler_attempted;
