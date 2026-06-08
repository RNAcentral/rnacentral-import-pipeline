LOAD CSV
FROM ALL FILENAMES MATCHING ~<qa.*csv$>
HAVING FIELDS (
  rna_id,
  upi,
  taxid,
  has_issue,
  incomplete_sequence,
  possible_contamination,
  missing_rfam_match,
  from_repetitive_region,
  possible_orf,
  possible_orf_stopfree,
  possible_orf_tcode,
  messages
)
INTO {{PGDATABASE}}?load_qa_status
TARGET COLUMNS (
  rna_id,
  upi,
  taxid,
  has_issue,
  incomplete_sequence,
  possible_contamination,
  missing_rfam_match,
  from_repetitive_region,
  possible_orf,
  possible_orf_stopfree,
  possible_orf_tcode,
  messages
)

WITH
    batch size = 32MB,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE qa_status
  ADD COLUMN IF NOT EXISTS possible_orf bool;
$$,
$$
ALTER TABLE qa_status
  ADD COLUMN IF NOT EXISTS possible_orf_stopfree bool;
$$,
$$
ALTER TABLE qa_status
  ADD COLUMN IF NOT EXISTS possible_orf_tcode bool;
$$,
$$
insert into qa_status (
  rna_id,
  upi,
  taxid,
  has_issue,
  incomplete_sequence,
  possible_contamination,
  missing_rfam_match,
  from_repetitive_region,
  possible_orf,
  possible_orf_stopfree,
  possible_orf_tcode,
  messages
) (
SELECT distinct on (rna_id)
  rna_id,
  upi,
  taxid,
  has_issue,
  incomplete_sequence,
  possible_contamination,
  missing_rfam_match,
  from_repetitive_region,
  possible_orf,
  possible_orf_stopfree,
  possible_orf_tcode,
  messages
FROM load_qa_status
)
ON CONFLICT (rna_id) DO UPDATE
SET
  has_issue = EXCLUDED.has_issue,
  incomplete_sequence = EXCLUDED.incomplete_sequence,
  possible_contamination = EXCLUDED.possible_contamination,
  missing_rfam_match = EXCLUDED.missing_rfam_match,
  from_repetitive_region = EXCLUDED.from_repetitive_region,
  possible_orf = EXCLUDED.possible_orf,
  possible_orf_stopfree = EXCLUDED.possible_orf_stopfree,
  possible_orf_tcode = EXCLUDED.possible_orf_tcode,
  messages = EXCLUDED.messages
;
$$,
$$
drop table load_qa_status;
$$
;
