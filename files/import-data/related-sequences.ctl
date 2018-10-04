LOAD CSV
FROM ALL FILENAMES MATCHING ~<related_sequences.*csv$>
HAVING FIELDS (
  source_accession,
  target_accession,
  relationship_type,
  methods
)
INTO {{PGDATABASE}}?load_rnc_related_sequences
TARGET COLUMNS (
  source_accession,
  target_accession,
  relationship_type,
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
  soruce_urs_taxid text,
  target_accession varchar(100) NOT NULL,
  relationship_type text NOT NULL,
  methods text[]
);
$$
;
