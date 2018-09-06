LOAD CSV
FROM ALL FILENAMES MATCHING ~<rfam_ontology_mappings.*csv$>
WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    ontology_term_id,
    rfam_model_id
)
INTO {{PGDATABASE}}?load_rfam_go_terms
TARGET COLUMNS
(
    ontology_term_id,
    rfam_model_id
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

AFTER LOAD DO
$$
insert into rfam_go_terms (
    ontology_term_id,
    rfam_model_id
) (
select
    ontology_term_id,
    rfam_model_id
from load_rfam_go_terms
)
ON CONFLICT (ontology_term_id, rfam_model_id) DO UPDATE SET
    ontology_term_id = excluded.ontology_term_id,
    rfam_model_id = excluded.rfam_model_id
;
$$,
$$
drop table load_rfam_go_terms;
$$
;
