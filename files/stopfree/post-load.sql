CREATE TABLE IF NOT EXISTS stopfree_results (
  urs_taxid TEXT PRIMARY KEY,
  stop_free_run_length integer,
  gc_content float,
  run_probability float,
  is_protein_coding bool
);

INSERT INTO stopfree_results (
  urs_taxid,
  stop_free_run_length,
  gc_content,
  run_probability,
  is_protein_coding
) (
SELECT
  urs_taxid,
  stop_free_run_length,
  gc_content,
  run_probability,
  is_protein_coding
FROM load_stopfree
)
ON CONFLICT (urs_taxid) DO UPDATE
SET
  stop_free_run_length = EXCLUDED.stop_free_run_length,
  gc_content = EXCLUDED.gc_content,
  run_probability = EXCLUDED.run_probability,
  is_protein_coding = EXCLUDED.is_protein_coding;
