LOAD CSV
FROM ALL FILENAMES MATCHING ~<cpat-orfs.*csv$>
HAVING FIELDS (
  urs,
  taxid,
  start_index,
  stop_index,
  metadata
)
INTO {{PGDATABASE}}?load_cpat_orfs
TARGET COLUMNS (
  urs,
  taxid,
  start_index,
  stop_index,
  metadata
)

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_cpat_orfs;
$$,
$$
CREATE TABLE load_cpat_orfs (
  urs TEXT NOT NULL,
  taxid int not null,
  start_index int not null,
  stop_index int not null,
  metadata jsonb not null
);
$$

AFTER LOAD DO
$$
DELETE FROM rnc_sequence_features features
USING load_cpat_orfs orfs
WHERE
  orfs.urs = features.upi
  and orfs.taxid = features.taxid
  and features.feature_name = 'cpat_orf'
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
  'cpat_orf',
  metadata
from load_cpat_orfs
);
$$,
$$
DROP TABLE load_cpat_orfs;
$$
;
