-- Minimal stand-in for the rnacen tables the SQL tests touch.
-- Only the columns those queries read are present; this is deliberately not a
-- copy of production DDL.

DROP TABLE IF EXISTS search_export_urs;
DROP TABLE IF EXISTS precompute_urs_taxid;
DROP TABLE IF EXISTS r2dt_results;
DROP TABLE IF EXISTS r2dt_models;
DROP TABLE IF EXISTS rna;

CREATE TABLE rna (
  urs text PRIMARY KEY,
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

-- Stand-ins for the tables files/import-data/post-release/001__regions.sql touches.

DROP MATERIALIZED VIEW IF EXISTS active_urs_taxids;
DROP MATERIALIZED VIEW IF EXISTS rnc_sequence_regions_active_provided;
DROP MATERIALIZED VIEW IF EXISTS rnc_sequence_regions_active_mapped;
DROP MATERIALIZED VIEW IF EXISTS rnc_sequence_regions_active;
DROP TABLE IF EXISTS load_rnc_sequence_regions;
DROP TABLE IF EXISTS xref;
DROP TABLE IF EXISTS ensembl_assembly;
DROP TABLE IF EXISTS rnc_gene_members;
DROP TABLE IF EXISTS rnc_gene_status;
DROP TABLE IF EXISTS rnc_genes;
DROP TABLE IF EXISTS rnc_accession_sequence_region;
DROP TABLE IF EXISTS rnc_sequence_exons;
DROP TABLE IF EXISTS rnc_sequence_regions;
DROP TABLE IF EXISTS rnc_rna_precomputed;

CREATE TABLE load_rnc_sequence_regions (
  accession text,
  urs_taxid text,
  region_name text NOT NULL,
  chromosome text,
  strand int,
  exon_start int,
  exon_stop int,
  assembly_id text,
  exon_count int
);

CREATE TABLE xref (
  ac text NOT NULL,
  urs text NOT NULL,
  taxid int NOT NULL,
  deleted char(1) NOT NULL
);

CREATE TABLE ensembl_assembly (
  assembly_id text PRIMARY KEY,
  assembly_ucsc text
);

CREATE TABLE rnc_sequence_regions (
  id serial PRIMARY KEY,
  urs_taxid text,
  region_name text NOT NULL,
  chromosome text,
  strand int,
  region_start int,
  region_stop int,
  assembly_id text,
  exon_count int,
  was_mapped bool,
  identity float
);
CREATE UNIQUE INDEX ON rnc_sequence_regions (md5(region_name), assembly_id);

CREATE TABLE rnc_sequence_exons (
  id serial PRIMARY KEY,
  region_id int NOT NULL REFERENCES rnc_sequence_regions(id),
  exon_start int,
  exon_stop int,
  UNIQUE (region_id, exon_start, exon_stop)
);

CREATE TABLE rnc_accession_sequence_region (
  accession text NOT NULL,
  region_id int NOT NULL REFERENCES rnc_sequence_regions(id),
  UNIQUE (accession, region_id)
);

CREATE TABLE rnc_genes (
  id serial PRIMARY KEY,
  member_count int,
  is_active bool
);

CREATE TABLE rnc_gene_members (
  rnc_gene_id int REFERENCES rnc_genes(id),
  locus_id int REFERENCES rnc_sequence_regions(id)
);

CREATE TABLE rnc_gene_status (
  region_id int REFERENCES rnc_sequence_regions(id)
);

CREATE TABLE rnc_rna_precomputed (
  urs_taxid text PRIMARY KEY,
  has_coordinates bool
);

CREATE MATERIALIZED VIEW active_urs_taxids AS
  SELECT DISTINCT urs || '_' || taxid AS urs_taxid FROM xref WHERE deleted = 'N';
CREATE MATERIALIZED VIEW rnc_sequence_regions_active_provided AS
  SELECT urs_taxid, assembly_id FROM rnc_sequence_regions WHERE NOT was_mapped;
CREATE MATERIALIZED VIEW rnc_sequence_regions_active_mapped AS
  SELECT urs_taxid, assembly_id FROM rnc_sequence_regions WHERE was_mapped;
CREATE MATERIALIZED VIEW rnc_sequence_regions_active AS
  SELECT urs_taxid, assembly_id FROM rnc_sequence_regions;
