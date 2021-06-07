LOAD CSV
FROM ALL FILENAMES MATCHING ~<precompute.*csv$>
HAVING FIELDS (
  id,
  upi,
  taxid [null if ""],
  is_active,
  description,
  rna_type,
  has_coordinates,
  databases,
  short_description,
  last_release,
  so_rna_type
)
INTO {{PGDATABASE}}?load_precomputed
TARGET COLUMNS (
  id,
  upi,
  taxid,
  is_active,
  description,
  rna_type,
  has_coordinates,
  databases,
  short_description,
  last_release,
  so_rna_type
)

WITH
    batch size = 32MB,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
insert into rnc_rna_precomputed (
  id,
  upi,
  taxid,
  is_active,
  rna_type,
  rfam_problems,
  description,
  databases,
  has_coordinates,
  short_description,
  last_release,
  so_rna_type
) (
SELECT DISTINCT ON (id)
  id,
  upi,
  taxid,
  is_active,
  rna_type,
  '{}'::jsonb,
  description,
  databases,
  has_coordinates,
  short_description,
  last_release,
  so_rna_type
FROM load_precomputed
)
ON CONFLICT (id) DO UPDATE
SET
  is_active = EXCLUDED.is_active,
  rna_type = EXCLUDED.rna_type,
  description = EXCLUDED.description,
  databases = EXCLUDED.databases,
  has_coordinates = EXCLUDED.has_coordinates,
  short_description = EXCLUDED.short_description,
  last_release = EXCLUDED.last_release,
  so_rna_type = EXCLUDED.so_rna_type
;
$$,
$$
drop table load_precomputed;
$$
;
