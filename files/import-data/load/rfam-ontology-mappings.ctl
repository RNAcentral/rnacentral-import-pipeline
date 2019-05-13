LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam_ontology_mappings.*csv$>
WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    rfam_model_id,
    ontology_term_id
)
INTO {{PGDATABASE}}?load_rfam_go_terms
TARGET COLUMNS
(
    rfam_model_id,
    ontology_term_id
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rfam_go_terms;
$$,
$$
create table load_rfam_go_terms (
    ontology_term_id character varying(10) COLLATE pg_catalog."default" NOT NULL,
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL
);
$$
;
