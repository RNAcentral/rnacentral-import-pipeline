LOAD CSV
FROM ALL FILENAMES MATCHING ~<pfam.*csv$>
HAVING FIELDS
(
  urs,
  sequence_start,
  sequence_stop,
  strand,
  pfam_model_id,
  model_start,
  model_stop,
  e_value,
  score
)
INTO {{PGDATABASE}}?load_pfam_model_hits
TARGET COLUMNS
(
  urs,
  sequence_start,
  sequence_stop,
  pfam_model_id,
  model_start,
  model_stop,
  e_value,
  score
)
WITH
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
SET work_mem='512MB';
$$,
$$
DROP TABLE IF EXISTS load_pfam_model_hits;
$$,
$$
CREATE TABLE load_pfam_model_hits (
  urs character varying(13) COLLATE pg_catalog."default" NOT NULL,
  pfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
  sequence_start integer NOT NULL,
  sequence_stop integer NOT NULL,
  sequence_completeness double precision,
  model_start integer NOT NULL,
  model_stop integer NOT NULL,
  model_completeness double precision,
  e_value double precision NOT NULL,
  score double precision NOT NULL
);
$$

AFTER LOAD DO
$$
CREATE INDEX ix__load_pfam_model_hits__urs ON load_pfam_model_hits (urs);
$$,
$$
CREATE INDEX ix__load_pfam_model_hits__pfam_model_id ON load_pfam_model_hits (pfam_model_id);
$$,
$$ INSERT INTO pfam_model_hits (
  urs,
  pfam_model_id,
  sequence_start,
  sequence_stop,
  sequence_completeness,
  model_start,
  model_stop,
  model_completeness,
  e_value,
  score
) (
SELECT
  load.urs,
  load.pfam_model_id,
  load.sequence_start,
  load.sequence_stop,
  abs((load.sequence_stop - load.sequence_start)::float) / rna.len::float,
  load.model_start,
  load.model_stop,
  abs((load.model_stop - load.model_start)::float) / models.length::float,
  load.e_value,
  load.score
FROM rna
JOIN pfam_models AS models,
JOIN load_pfam_model_hits AS load
WHERE
  rna.upi = load.urs
  AND models.pfam_model_id = load.pfam_model_id
)
ON CONFLICT (sequence_start, sequence_stop, model_start, model_stop, pfam_model_id, urs)
DO UPDATE
SET
  e_value = excluded.e_value,
  score = excluded.score
;
$$,
$$
INSERT INTO pfam_analyzed_sequences (
    urs,
    date,
    total_matches,
    total_family_matches
) (
SELECT
    urs,
    current_date date,
    count(*) total_matches,
    count(distinct pfam_model_id) total_family_matches
FROM pfam_model_hits
GROUP BY urs
)
ON CONFLICT (urs)
DO UPDATE
SET
  total_matches = EXCLUDED.total_matches,
  total_family_matches = EXCLUDED.total_family_matches
;
$$,
$$
drop table load_pfam_model_hits;
$$
;
