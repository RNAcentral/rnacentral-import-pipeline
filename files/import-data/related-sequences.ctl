LOAD CSV
FROM stdin
HAVING FIELDS (
    source_accession,
    target_accession,
    relationship,
    methods
)
INTO {{PGDATABASE}}?load_rnc_related_sequences
TARGET COLUMNS (
    source_accession,
    target_accession,
    relationship,
    methods
)

WITH truncate,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rnc_related_sequences;
$$,
$$
create table load_rnc_related_sequences (
    source_accession varchar(100) NOT NULL,
    target_accession varchar(100) NOT NULL,
    relationship text NOT NULL,
    methods text[]
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_related_sequences (
    source_accession,
    target_accession,
    relationship,
    methods,
    is_internal
) (
select
    load.source_accession,
    load.target_accession,
    load.relationship::related_sequence_relationship,
    load.methods,
    exists (select 1 from rnc_accessions where accession = load.target_accession)
from load_rnc_related_sequences load
)
ON CONFLICT (source_accession, target_accession, relationship) DO UPDATE
SET
    methods = EXCLUDED.methods || rnc_related_sequences.methods,
    is_internal = EXCLUDED.is_internal or rnc_related_sequences.is_internal
;
$$,
$$
drop table load_rnc_related_sequences;
$$
;
