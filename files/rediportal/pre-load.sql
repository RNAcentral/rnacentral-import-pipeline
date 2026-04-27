-- Ensure the rediportal staging table exists before bin/load-parquet writes
-- to it. Replaces the BEFORE LOAD DO block in files/rediportal/load.ctl,
-- which in the CSV path drops+recreates the table every run. With parquet,
-- bin/load-parquet --truncate handles the wipe; this script just guarantees
-- the table is there. Idempotent: safe to run on every pipeline invocation.
--
-- Mirrors the column types from the original ctl: nullable on every column,
-- jsonb metadata. The parser populates everything except ``accession``.

CREATE TABLE IF NOT EXISTS load_rediportal_features (
    upi text,
    taxid int,
    accession text,
    start int,
    stop int,
    feature_name text,
    metadata jsonb,
    feature_provider text
);
