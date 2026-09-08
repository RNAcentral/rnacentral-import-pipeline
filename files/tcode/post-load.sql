CREATE TABLE IF NOT EXISTS tcode_results (
  urs_taxid TEXT PRIMARY KEY,
  length integer,
  mean_score float,
  std_score float,
  is_protein_coding bool
);

INSERT INTO tcode_results (
  urs_taxid,
  length,
  mean_score,
  std_score,
  is_protein_coding
) (
SELECT
  urs_taxid,
  length,
  mean_score,
  std_score,
  is_protein_coding
FROM load_tcode
)
ON CONFLICT (urs_taxid) DO UPDATE
SET
  length = EXCLUDED.length,
  mean_score = EXCLUDED.mean_score,
  std_score = EXCLUDED.std_score,
  is_protein_coding = EXCLUDED.is_protein_coding;
