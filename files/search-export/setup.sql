CREATE TABLE search_export_urs (
  id BIGSERIAL PRIMARY KEY,
  urs TEXT NOT NULL,
  taxid INT NOT NULL,
  urs_taxid TEXT NOT NULL
);

CREATE TABLE search_export_accessions (
  id BIGSERIAL PRIMARY KEY,
  search_export_id BIGSERIAL NOT NULL,
  urs_taxid TEXT NOT NULL,
  accession TEXT NOT NULL,
  classification TEXT,
  common_name TEXT,
  database TEXT,
  external_id TEXT,
  function TEXT,
  gene TEXT,
  gene_synonym TEXT,
  genes TEXT,
  locus_tag TEXT,
  locus_tags TEXT,
  non_coding_id TEXT,
  note TEXT,
  optional_id TEXT,
  organelle TEXT,
  organelles TEXT,
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
) ON CONFLICT DO NOTHING;

CREATE INDEX ix__search_todo__urs_taxid ON search_export_urs(urs, taxid);
