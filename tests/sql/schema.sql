-- Minimal stand-in for the rnacen tables the r2dt display queries touch.
-- Only the columns those queries read are present; this is deliberately not a
-- copy of production DDL.

DROP TABLE IF EXISTS search_export_urs;
DROP TABLE IF EXISTS precompute_urs_taxid;
DROP TABLE IF EXISTS r2dt_results;
DROP TABLE IF EXISTS r2dt_models;
DROP TABLE IF EXISTS rna;

CREATE TABLE rna (
  upi text PRIMARY KEY,
  len int NOT NULL
);

CREATE TABLE r2dt_models (
  id int PRIMARY KEY,
  model_name text NOT NULL,
  model_source text NOT NULL,
  model_length int,
  model_basepair_count int,
  so_term_id text
);

CREATE TABLE r2dt_results (
  urs text NOT NULL,
  model_id int NOT NULL REFERENCES r2dt_models(id),
  inferred_should_show bool,
  assigned_should_show bool,
  sequence_coverage float,
  model_coverage float,
  basepair_count int
);

CREATE TABLE search_export_urs (
  id bigint PRIMARY KEY,
  urs text NOT NULL,
  urs_taxid text NOT NULL
);

CREATE TABLE precompute_urs_taxid (
  id bigint PRIMARY KEY,
  precompute_urs_id bigint NOT NULL,
  urs text NOT NULL,
  urs_taxid text NOT NULL
);
