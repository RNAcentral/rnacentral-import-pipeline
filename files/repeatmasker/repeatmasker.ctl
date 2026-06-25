LOAD CSV
FROM ALL FILENAMES MATCHING ~<repeatmasker-results.*csv$>
HAVING FIELDS (
  urs_taxid,
  has_repeats,
  repeat_coverage,
  repeat_count
)
INTO {{PGDATABASE}}?load_repeatmasker
TARGET COLUMNS (
  urs_taxid,
  has_repeats,
  repeat_coverage,
  repeat_count
)

BEFORE LOAD DO
$$
CREATE TABLE IF NOT EXISTS repeatmasker_results (
  urs_taxid TEXT PRIMARY KEY,
  has_repeats bool,
  repeat_coverage float,
  repeat_count integer
);
$$

AFTER LOAD DO
$$
INSERT INTO repeatmasker_results (
  urs_taxid,
  has_repeats,
  repeat_coverage,
  repeat_count
) (
SELECT
  urs_taxid,
  has_repeats,
  repeat_coverage,
  repeat_count
from load_repeatmasker
) ON CONFLICT (urs_taxid) DO UPDATE
SET
  has_repeats = EXCLUDED.has_repeats,
  repeat_coverage = EXCLUDED.repeat_coverage,
  repeat_count = EXCLUDED.repeat_count
;
$$,
$$
DROP TABLE load_repeatmasker;
$$
;
