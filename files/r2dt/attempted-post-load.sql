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
