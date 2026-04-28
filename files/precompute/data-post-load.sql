INSERT INTO rnc_rna_precomputed (
  id,
  upi,
  taxid,
  is_active,
  rna_type,
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
  so_rna_type = EXCLUDED.so_rna_type;

DROP TABLE load_precomputed;
