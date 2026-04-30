LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam.*csv$>
HAVING FIELDS
(
  UPI,
  SEQUENCE_START,
  SEQUENCE_STOP,
  STRAND,
  RFAM_MODEL_ID,
  MODEL_START,
  MODEL_STOP,
  OVERLAP,
  E_VALUE,
  SCORE
)
INTO {{PGDATABASE}}?load_rfam_model_hits
TARGET COLUMNS
(
  UPI,
  SEQUENCE_START, SEQUENCE_STOP,
  RFAM_MODEL_ID,
  MODEL_START,
  MODEL_STOP,
  OVERLAP,
  E_VALUE,
  SCORE
)
WITH fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
set work_mem='512MB';
$$

AFTER LOAD DO
$$
create index ix__load_rfam_model_hits__upi on load_rfam_model_hits (upi);
$$,
$$
create index ix__load_rfam_model_hits__rfam_model_id on load_rfam_model_hits (rfam_model_id);
$$,
$$
DELETE FROM rfam_model_hits hits
USING load_rfam_model_hits load
WHERE
  hits.upi = load.upi
;
$$,
$$
insert into rfam_model_hits (
  sequence_start,
  sequence_stop,
  sequence_completeness,
  model_start,
  model_stop,
  model_completeness,
  overlap,
  e_value,
  score,
  rfam_model_id,
  upi
) (
select
  load.sequence_start,
  load.sequence_stop,
  abs((load.sequence_stop - load.sequence_start)::float) / rna.len::float,
  load.model_start,
  load.model_stop,
  abs((load.model_stop - load.model_start)::float) / models.length::float,
  load.overlap,
  load.e_value,
  load.score,
  load.rfam_model_id,
  load.upi
from rna, rfam_models as models, load_rfam_model_hits as load
where
  rna.upi = load.upi
  and models.rfam_model_id = load.rfam_model_id
)
ON CONFLICT (sequence_start, sequence_stop, model_start, model_stop, rfam_model_id, upi)
DO UPDATE
SET
  e_value = excluded.e_value,
  score = excluded.score,
  overlap = excluded.overlap
;
$$,
$$
drop index ix__load_rfam_model_hits__rfam_model_id;
$$,
$$
drop index ix__load_rfam_model_hits__upi;
$$,
$$
drop table load_rfam_model_hits;
$$
;
