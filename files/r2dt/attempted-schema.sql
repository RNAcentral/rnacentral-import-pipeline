-- Staging table for the r2dt "attempted" side.
--
-- Kept in step with the BEFORE LOAD DO body of files/r2dt/attempted.ctl: the
-- pgloader path runs that, the parquet path runs this. bin/load-parquet does
-- not execute a ctl's BEFORE LOAD DO, so without this file the parquet branch
-- depends on files/schema/create_load.sql having been run by a previous
-- import-data release -- a silent cross-workflow dependency.

DROP TABLE IF EXISTS load_traveler_attempted;

CREATE UNLOGGED TABLE load_traveler_attempted (
  urs text primary key,
  r2dt_version text
);
