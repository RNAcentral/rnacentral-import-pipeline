LOAD CSV
FROM 'term-info.csv'
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

AFTER LOAD DO
$$ insert into ontology_terms (
  ontology_term_id,
  ontology,
  name,
  definition
) (
select distinct
  ontology_term_id,
  ontology,
  name,
  definition
from load_ontology_terms
)
ON CONFLICT (ontology_term_id) DO UPDATE SET
  ontology_term_id = excluded.ontology_term_id,
  ontology = excluded.ontology,
  name = excluded.name,
  definition = excluded.definition
;
$$,
$$
drop table load_ontology_terms;
$$
;
