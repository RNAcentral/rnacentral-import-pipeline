LOAD CSV
FROM ALL FILENAMES MATCHING ~<attempted.*csv$>
HAVING FIELDS
(
  urs,
  model_source,
  source_version
)
INTO {{PGDATABASE}}?load_qa_rfam_attempted
TARGET COLUMNS
(
  urs,
  model_source,
  source_version
)
WITH fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
set work_mem='512MB';
$$

AFTER LOAD DO
$$
INSERT INTO pipeline_tracking_qa_rfam (
  urs,
  model_source,
  source_version,
  last_run
) (
select 
  urs,
  model_source,
  source_version,
  last_run
FROM load_qa_rfam_attempted
) ON CONFLICT (urs, model_source) DO UPDATE
SET 
  model_source = EXCLUDED.model_source,
  source_version = EXCLUDED.source_version,
  last_run = EXCLUDED.last_run
$$
;

