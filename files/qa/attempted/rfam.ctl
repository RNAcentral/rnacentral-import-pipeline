LOAD CSV
FROM ALL FILENAMES MATCHING ~<attempted.*csv$>
HAVING FIELDS
(
  urs,
  model_source,
  source_version,
  last_done
)
INTO {{PGDATABASE}}?load_qa_rfam_attempted
TARGET COLUMNS
(
  urs,
  model_source,
  source_version,
  last_done
)
WITH fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
set work_mem='512MB';
$$,
$$

AFTER LOAD DO
$$
INSERT INTO pipeline_tracking_qa_rfam (
  urs,
  model_source,
  source_version,
  last_done
) (
select 
  urs,
  model_source,
  source_version,
  last_done
FROM load_qa_rfam_attempted
) ON CONFLICT (urs) DO UPDATE
SET 
  qa_analysis = EXCLUDED.qa_analysis,
  last_done = EXCLUDED.last_done
$$
;

