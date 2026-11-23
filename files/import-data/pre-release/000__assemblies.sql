BEGIN;

-- Ensembl dropped its subdomains, so new rows leave this blank.
ALTER TABLE ensembl_assembly ALTER COLUMN subdomain DROP NOT NULL;

INSERT INTO ensembl_assembly (
  assembly_id,
  assembly_full_name,
  gca_accession,
  assembly_ucsc,
  common_name,
  taxid,
  ensembl_url,
  blat_mapping,
  example_chromosome,
  example_end,
  example_start
) (
SELECT DISTINCT ON (assembly_id)
  assembly_id,
  assembly_full_name,
  gca_accession,
  assembly_ucsc,
  common_name,
  taxid,
  ensembl_url,
  blat_mapping,
  example_chromosome,
  example_end,
  example_start
FROM load_assemblies
) ON CONFLICT (assembly_id) DO UPDATE
-- ensembl_url is never updated: species.json renamed species we already hold.
-- Nor is taxid: species.json can give a strain where the files give the species.
SET
  assembly_full_name = EXCLUDED.assembly_full_name,
  gca_accession = EXCLUDED.gca_accession,
  assembly_ucsc = COALESCE(EXCLUDED.assembly_ucsc, ensembl_assembly.assembly_ucsc),
  common_name = EXCLUDED.common_name,
  blat_mapping = EXCLUDED.blat_mapping,
  example_chromosome = EXCLUDED.example_chromosome,
  example_end = EXCLUDED.example_end,
  example_start = EXCLUDED.example_start
;

DROP TABLE load_assemblies;

COMMIT;
