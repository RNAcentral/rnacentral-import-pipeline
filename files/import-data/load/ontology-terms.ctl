LOAD CSV
FROM 'ontology_terms.csv'
WITH ENCODING ISO-8859-14
HAVING FIELDS
(
  ontology_term_id,
  ontology,
  name,
  definition
)
INTO {{PGDATABASE}}?load_ontology_terms
TARGET COLUMNS
(
  ontology_term_id,
  ontology,
  name,
  definition
)

WITH
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_ontology_terms;
$$,
$$
create table load_ontology_terms (
  ontology_term_id varchar(15),
  ontology varchar(5),
  name text,
  definition text
);
$$
;
