LOAD CSV
FROM ALL FILENAMES MATCHING ~<cpat-results.*csv$>
HAVING FIELDS (
  urs_taxid,
  fickett_score,
  hexamer_score,
  coding_probability,
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

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_cpat;
$$,
$$
CREATE TABLE load_cpat (
  urs_taxid TEXT not null,
  fickett_score float not null,
  hexamer_score float not null,
  coding_probability float not null,
  is_protein_coding bool not null
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_cpat_results (
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
