INSERT INTO cpat_results (
  urs_taxid,
  fickett_score,
  hexamer_score,
  coding_probability,
  is_protein_coding
) (
SELECT
  urs_taxid,
  fickett_score,
  hexamer_score,
  coding_probability,
  is_protein_coding
FROM load_cpat
)
ON CONFLICT (urs_taxid) DO UPDATE
SET
  fickett_score = EXCLUDED.fickett_score,
  hexamer_score = EXCLUDED.hexamer_score,
  coding_probability = EXCLUDED.coding_probability,
  is_protein_coding = EXCLUDED.is_protein_coding;
