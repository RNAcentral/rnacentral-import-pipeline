DELETE FROM rnc_sequence_features features
USING load_cpat_orfs orfs
WHERE
  orfs.urs = features.urs
  AND orfs.taxid = features.taxid
  AND features.feature_name = 'cpat_orf';

INSERT INTO rnc_sequence_features (
  urs,
  taxid,
  "start",
  "stop",
  feature_name,
  metadata
) (
SELECT
  urs,
  taxid,
  start_index,
  stop_index,
  'cpat_orf',
  metadata
FROM load_cpat_orfs
)
ON CONFLICT (urs, taxid, accession, start, stop, feature_name) DO UPDATE
SET
  metadata = excluded.metadata;
