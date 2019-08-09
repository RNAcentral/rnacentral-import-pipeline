LOAD CSV
FROM ALL FILENAMES MATCHING ~<go_publication_mappings.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
    rna_id,
    qualifier,
    ontology_term_id,
    assigned_by,
    evidence_code,
    pubmed_id
)
INTO {{PGDATABASE}}?load_go_term_publication_map
TARGET COLUMNS (
    rna_id,
    qualifier,
    ontology_term_id,
    assigned_by,
    evidence_code,
    pubmed_id
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_go_term_publication_map;
$$,
$$
create table if not exists load_go_term_publication_map (
    rna_id varchar(50),
    qualifier text,
    assigned_by varchar(50),
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    pubmed_id text
);
$$
AFTER LOAD DO
$$
update load_go_term_publication_map
set pubmed_id = trim(leading 'pmid:' from pubmed_id)
$$
;
