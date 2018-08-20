LOAD CSV
FROM annotations.csv
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
    qualifier varchar(30),
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    assigned_by varchar(50),
    extensions jsonb
);
$$

AFTER LOAD DO
$$
INSERT INTO go_term_annotations (
    rna_id,
    qualifier,
    ontology_term_id,
    evidence_code,
    assigned_by,
    extensions
) (
SELECT
    rna_id,
    qualifier,
    ontology_term_id,
    evidence_code,
    assigned_by,
    extensions
FROM load_go_term_annotations
)
ON CONFLICT (rna_id, qualifier, ontology_term_id, evidence_code, assigned_by)
DO UPDATE
SET
    extensions = excluded.extensions
;
$$,
$$
DROP TABLE load_go_term_annotations;
$$
;
