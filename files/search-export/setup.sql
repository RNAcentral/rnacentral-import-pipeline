BEGIN TRANSACTION;

DROP TABLE IF EXISTS search_export_urs;
CREATE TABLE IF NOT EXISTS search_export_urs (
  id BIGSERIAL PRIMARY KEY,
  urs TEXT NOT NULL,
  taxid INT NOT NULL,
  urs_taxid TEXT NOT NULL
);

DROP TABLE IF EXISTS search_export_accessions;
CREATE TABLE IF NOT EXISTS search_export_accessions (
  id BIGSERIAL PRIMARY KEY,
  search_export_id BIGSERIAL NOT NULL,
  urs TEXT NOT NULL,
  taxid INT NOT NULL,
  urs_taxid TEXT NOT NULL,
  accession TEXT NOT NULL,
  lineage TEXT,
  common_name TEXT,
  database TEXT,
  external_id TEXT,
  function TEXT,
  gene TEXT,
  gene_synonym TEXT,
  locus_tag TEXT,
  non_coding_id TEXT,
  note TEXT,
  optional_id TEXT,
  organelle TEXT,
  parent_accession TEXT,
  product TEXT,
  species TEXT,
  standard_name TEXT
);

INSERT INTO search_export_urs (urs_taxid, urs, taxid) (
  SELECT
    pre.id,
    pre.upi,
    pre.taxid
  FROM rnc_rna_precomputed pre
  WHERE
    pre.is_active = true
    and taxid is not null
) ON CONFLICT DO NOTHING;

CREATE INDEX ix__search_todo__urs_taxid ON search_export_urs(urs, taxid);

COMMIT;
