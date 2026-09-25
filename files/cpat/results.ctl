LOAD CSV
FROM ALL FILENAMES MATCHING ~<cpat-results.*csv$>
HAVING FIELDS (
  urs_taxid,
  fickett_score [null if ""],
  hexamer_score [null if ""],
  coding_probability [null if ""],
  is_protein_coding
)
INTO {{PGDATABASE}}?load_cpat
TARGET COLUMNS (
  urs_taxid,
  fickett_score,
  hexamer_score,
  coding_probability,
  is_protein_coding
)

-- Sequences without an ORF are stored with null scores; the tables predate that.
BEFORE LOAD DO
$$
ALTER TABLE load_cpat
  ALTER COLUMN fickett_score DROP NOT NULL,
  ALTER COLUMN hexamer_score DROP NOT NULL,
  ALTER COLUMN coding_probability DROP NOT NULL;
$$,
$$
ALTER TABLE cpat_results
  ALTER COLUMN fickett_score DROP NOT NULL,
  ALTER COLUMN hexamer_score DROP NOT NULL,
  ALTER COLUMN coding_probability DROP NOT NULL;
$$

AFTER LOAD DO
$$
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
from load_cpat
) ON CONFLICT (urs_taxid) DO UPDATE
SET
  fickett_score = EXCLUDED.fickett_score,
  hexamer_score = EXCLUDED.hexamer_score,
  coding_probability = EXCLUDED.coding_probability,
  is_protein_coding = EXCLUDED.is_protein_coding
;
$$,
$$
DROP TABLE load_cpat;
$$
;
