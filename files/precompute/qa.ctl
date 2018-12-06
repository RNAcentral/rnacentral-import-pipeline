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
  messages
)

WITH
    batch size = 32MB,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_qa_status (
  rna_id varchar(44) NOT NULL,
  upi varchar(26) NOT NULL,
  taxid int8 NOT NULL,
  has_issue bool,
  incomplete_sequence bool,
  possible_contamination bool,
  missing_rfam_match bool,
  messages jsonb
);
$$,
$$
truncate table load_qa_status;
$$

AFTER LOAD DO
$$
insert into qa_status (
  rna_id,
  upi,
  taxid,
  has_issue,
  incomplete_sequence,
  possible_contamination,
  missing_rfam_match,
  messages
) (
SELECT distinct
  rna_id,
  upi,
  taxid,
  has_issue,
  incomplete_sequence,
  possible_contamination,
  missing_rfam_match,
  messages
FROM load_qa_status
)
ON CONFLICT (rna_id) DO UPDATE
SET
  has_issue = EXCLUDED.has_issue,
  incomplete_sequence = EXCLUDED.incomplete_sequence,
  possible_contamination = EXCLUDED.possible_contamination,
  missing_rfam_match = EXCLUDED.missing_rfam_match,
  messages = EXCLUDED.messages
;
$$,
$$
drop table load_qa_status;
$$
;
