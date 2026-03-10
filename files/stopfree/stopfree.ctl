LOAD CSV
FROM ALL FILENAMES MATCHING ~<stopfree-results.*csv$>
HAVING FIELDS (
  urs_taxid,
  length,
  stop_free_run_length,
  gc_content,
  run_probability,
  is_protein_coding
)
INTO {{PGDATABASE}}?load_stopfree
TARGET COLUMNS (
  urs_taxid,
  length,
  stop_free_run_length,
  gc_content,
  run_probability,
  is_protein_coding
)

BEFORE LOAD DO
$$
CREATE TABLE IF NOT EXISTS stopfree_results (
  urs_taxid TEXT PRIMARY KEY,
  length integer,
  stop_free_run_length integer,
  gc_content float,
  run_probability float,
  is_protein_coding bool
);
$$,
$$
DROP TABLE IF EXISTS load_stopfree;
$$,
$$
CREATE TABLE load_stopfree (
  urs_taxid TEXT not null,
  length integer,
  stop_free_run_length integer,
  gc_content float,
  run_probability float,
  is_protein_coding bool
);
$$

AFTER LOAD DO
$$
INSERT INTO stopfree_results (
  urs_taxid,
  length,
  stop_free_run_length,
  gc_content,
  run_probability,
  is_protein_coding
) (
SELECT
  urs_taxid,
  length,
  stop_free_run_length,
  gc_content,
  run_probability,
  is_protein_coding
from load_stopfree
) ON CONFLICT (urs_taxid) DO UPDATE
SET
  length = EXCLUDED.length,
  stop_free_run_length = EXCLUDED.stop_free_run_length,
  gc_content = EXCLUDED.gc_content,
  run_probability = EXCLUDED.run_probability,
  is_protein_coding = EXCLUDED.is_protein_coding
;
$$,
$$
DROP TABLE load_stopfree;
$$
;
