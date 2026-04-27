-- Ensure the rfam-scan staging table exists before bin/load-parquet writes
-- to it. Replaces the BEFORE LOAD DO block in files/rfam-scan/load.ctl, which
-- in the CSV path drops+recreates the table every run. With parquet,
-- bin/load-parquet --truncate handles the wipe; this script just guarantees
-- the table is there. Idempotent: safe to run on every pipeline invocation.

CREATE TABLE IF NOT EXISTS load_rfam_model_hits (
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
