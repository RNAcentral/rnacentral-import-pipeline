LOAD CSV
FROM ALL FILENAMES MATCHING ~<tcode-results.*csv$>
HAVING FIELDS (
  urs_taxid,
  length,
  mean_score,
  std_score,
  is_protein_coding
)
INTO {{PGDATABASE}}?load_tcode
TARGET COLUMNS (
  urs_taxid,
  length,
  mean_score,
  std_score,
  is_protein_coding
)

BEFORE LOAD DO
$$
CREATE TABLE IF NOT EXISTS tcode_results (
  urs_taxid TEXT PRIMARY KEY,
  length integer,
  mean_score float,
  std_score float,
  is_protein_coding bool
);
$$,
$$
DROP TABLE IF EXISTS load_tcode;
$$,
$$
CREATE TABLE load_tcode (
  urs_taxid TEXT not null,
  length integer,
  mean_score float,
  std_score float,
  is_protein_coding bool
);
$$

AFTER LOAD DO
$$
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
from load_tcode
) ON CONFLICT (urs_taxid) DO UPDATE
SET
  length = EXCLUDED.length,
  mean_score = EXCLUDED.mean_score,
  std_score = EXCLUDED.std_score,
  is_protein_coding = EXCLUDED.is_protein_coding
;
$$,
$$
DROP TABLE load_tcode;
$$
;
