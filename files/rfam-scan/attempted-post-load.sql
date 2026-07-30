SET work_mem='512MB';

INSERT INTO pipeline_tracking_qa_scan (
  urs,
  model_source,
  source_version,
  last_run
) (
SELECT DISTINCT ON (urs)
  urs,
  model_source,
  source_version,
  last_run
FROM load_qa_rfam_attempted
) ON CONFLICT (urs, model_source) DO UPDATE
SET
  model_source = EXCLUDED.model_source,
  source_version = EXCLUDED.source_version,
  last_run = EXCLUDED.last_run;
