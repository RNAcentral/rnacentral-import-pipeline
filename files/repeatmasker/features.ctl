LOAD CSV
FROM ALL FILENAMES MATCHING ~<repeatmasker-features.*csv$>
HAVING FIELDS (
  urs,
  taxid,
  start_index,
  stop_index,
  metadata
)
INTO {{PGDATABASE}}?load_repeatmasker_features
TARGET COLUMNS (
  urs,
  taxid,
  start_index,
  stop_index,
  metadata
)

AFTER LOAD DO
$$
DELETE FROM rnc_sequence_features features
USING load_repeatmasker_features regions
WHERE
  regions.urs = features.upi
  and regions.taxid = features.taxid
  and features.feature_name = 'repeatmasker_region'
;
$$,
$$
INSERT INTO rnc_sequence_features (
  upi,
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
  'repeatmasker_region',
  metadata
from load_repeatmasker_features
) ON CONFLICT (upi, taxid, accession, start, stop, feature_name) DO UPDATE
SET
  metadata = excluded.metadata;
$$,
$$
DROP TABLE load_repeatmasker_features;
$$
;
