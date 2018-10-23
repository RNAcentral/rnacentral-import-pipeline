LOAD CSV
FROM ALL FILENAMES MATCHING ~<publications.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
    ref_pubmed_id,
    authors,
    location,
    title,
    doi
)
INTO {{PGDATABASE}}?load_ref_pubmed
TARGET COLUMNS (
    ref_pubmed_id,
    authors,
    location,
    title,
    doi
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_ref_pubmed;
$$,
$$
create table load_ref_pubmed (
    ref_pubmed_id int,
    authors text,
    location text,
    title text,
    doi text
);
$$
;
