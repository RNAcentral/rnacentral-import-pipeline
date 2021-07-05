LOAD CSV
FROM ALL FILENAMES MATCHING ~<cpat-orfs.*csv$>
HAVING FIELDS (
  urs_taxid,
  start_index,
  stop_index,
  strand
)
INTO {{PGDATABASE}}?load_cpat_orfs
TARGET COLUMNS (
  urs_taxid,
  start_index,
  stop_index,
  strand
)

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS load_cpat_orfs;
$$,
$$
CREATE TABLE load_cpat_orfs (
  urs_taxid TEXT not null,
  start_index int not null,
  stop_index int not null,
  strand strand not null
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_cpat_orfs (
  urs_taxid,
  start_index,
  stop_index,
  strand
) (
SELECT
  urs_taxid,
  start_index,
  stop_index,
  strand
from load_cpat_orfs
) ON CONFLICT (urs_taxid) DO UPDATE
SET
  start_index = EXCLUDED.start_index,
  stop_index = EXCLUDED.stop_index,
  strand = EXCLUDED.strand
;
$$,
$$
DROP TABLE load_cpat_orfs;
$$
;

