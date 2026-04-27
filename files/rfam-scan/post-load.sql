-- Merge fresh hits from load_rfam_model_hits into rfam_model_hits.
-- Lifted from the AFTER LOAD DO body of files/rfam-scan/load.ctl. Run by
-- bin/load-parquet --post-load after the staging table is populated; the
-- entire script executes in a single transaction so a failure leaves
-- rfam_model_hits untouched.
--
-- The completeness columns are computed here (not at load time) because they
-- depend on joins against rna and rfam_models — data that's only present in
-- the production database, not in the parquet input.

CREATE INDEX IF NOT EXISTS ix__load_rfam_model_hits__upi
  ON load_rfam_model_hits (upi);
CREATE INDEX IF NOT EXISTS ix__load_rfam_model_hits__rfam_model_id
  ON load_rfam_model_hits (rfam_model_id);

DELETE FROM rfam_model_hits hits
USING load_rfam_model_hits load
WHERE hits.upi = load.upi;

INSERT INTO rfam_model_hits (
  sequence_start,
  sequence_stop,
  sequence_completeness,
  model_start,
  model_stop,
  model_completeness,
  overlap,
  e_value,
  score,
  rfam_model_id,
  upi
)
SELECT
  load.sequence_start,
  load.sequence_stop,
  abs((load.sequence_stop - load.sequence_start)::float) / rna.len::float,
  load.model_start,
  load.model_stop,
  abs((load.model_stop - load.model_start)::float) / models.length::float,
  load.overlap,
  load.e_value,
  load.score,
  load.rfam_model_id,
  load.upi
FROM rna, rfam_models AS models, load_rfam_model_hits AS load
WHERE rna.upi = load.upi
  AND models.rfam_model_id = load.rfam_model_id
ON CONFLICT (sequence_start, sequence_stop, model_start, model_stop, rfam_model_id, upi)
DO UPDATE SET
  e_value = excluded.e_value,
  score = excluded.score,
  overlap = excluded.overlap;

DROP INDEX IF EXISTS ix__load_rfam_model_hits__rfam_model_id;
DROP INDEX IF EXISTS ix__load_rfam_model_hits__upi;
