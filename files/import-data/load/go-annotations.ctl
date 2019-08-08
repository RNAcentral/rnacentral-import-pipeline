LOAD CSV
FROM ALL FILENAMES MATCHING ~<annotations.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
  rna_id,
  qualifier,
  ontology_term_id,
  evidence_code,
  extensions,
  assigned_by
)
INTO {{PGDATABASE}}?load_go_term_annotations
TARGET COLUMNS (
  rna_id,
  qualifier,
  ontology_term_id,
  evidence_code,
  extensions,
  assigned_by
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
DROP TABLE if exists load_go_term_annotations;
$$,
$$
create table load_go_term_annotations (
    rna_id varchar(50),
    qualifier,
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    assigned_by varchar(50),
    extensions jsonb
);
$$
;
