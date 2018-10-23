INSERT INTO protein_info (
  protein_accession,
  description,
  label,
  synonyms
) (
SELECT DISTINCT
  protein_accession,
  description,
  label,
  synonyms
FROM load_protein_info
)
ON CONFLICT (protein_accession) DO UPDATE
SET
  description = EXCLUDED.description,
  label = EXCLUDED.label,
  synonyms = EXCLUDED.synonyms
;

DROP TABLE load_protein_info;
