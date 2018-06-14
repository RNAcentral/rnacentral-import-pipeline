LOAD CSV
FROM 'hits.csv'
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
INTO {{PGDATABASE}}
TARGET COLUMNS
(
  UPI,
  SEQUENCE_START,
  SEQUENCE_STOP,
  RFAM_MODEL_ID,
  MODEL_START,
  MODEL_STOP,
  OVERLAP,
  E_VALUE,
  SCORE
)
WITH skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
  set work_mem='512MB';
$$,
$$
create table if not exists load_rfam_model_hits (
    sequence_start integer NOT NULL,
    sequence_stop integer NOT NULL,
    sequence_completeness double precision,
    model_start integer NOT NULL,
    model_stop integer NOT NULL,
    model_completeness double precision,
    overlap character varying(30) COLLATE pg_catalog."default" NOT NULL,
    e_value double precision NOT NULL,
    score double precision NOT NULL,
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    upi character varying(13) COLLATE pg_catalog."default" NOT NULL
);
$$,
$$
truncate table load_rfam_model_hits;
$$

AFTER LOAD DO
$$
create index ix__load_rfam_model_hits__upi on load_rfam_model_hits (upi);
$$,
$$
create index ix__load_rfam_model_hits__rfam_model_id on load_rfam_model_hits (rfam_model_id);
$$,
$$ insert into rfam_model_hits (
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
);
$$,
$$
update rfam_analyzed_sequences
set
    total_matches = counts.total,
    total_family_matches = counts.family_count
from (
    select
        upi,
        count(*) total,
        count(distinct rfam_model_id) family_count
    from rfam_model_hits
    group by upi
) as counts
where
    rfam_analyzed_sequences.upi = counts.upi
;
$$,
$$
drop index ix__load_rfam_model_hits__rfam_model_id;
$$,
$$
drop index ix__load_rfam_model_hits__upi;
$$,
$$
truncate table load_rfam_model_hits;
$$,
$$
drop table load_rfam_model_hits;
$$
;
