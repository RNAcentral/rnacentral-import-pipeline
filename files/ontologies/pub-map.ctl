LOAD CSV
FROM ALL FILENAMES MATCHING ~<pub_map.*csv$>
WITH ENCODING ISO-8859-14

HAVING FIELDS (
    rna_id,
    qualifier,
    ontology_term_id,
    assigned_by,
    evidence_code,
    pmid
)
INTO {{PGDATABASE}}?load_go_term_publication_map
TARGET COLUMNS (
    rna_id,
    qualifier,
    ontology_term_id,
    assigned_by,
    evidence_code,
    pmid
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
    qualifier varchar(30),
    assigned_by varchar(50),
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    pubmed_id int
);
$$

AFTER LOAD DO
$$
insert into go_term_publication_map (
    go_term_annotation_id,
    ref_pubmed_id
) (
select distinct
    annotations.go_term_annotation_id,
    load_go_term_publication_map.pubmed_id
from load_go_term_publication_map
join go_term_annotations annotations
on
    annotations.rna_id = load_go_term_publication_map.rna_id
    AND annotations.qualifier = load_go_term_publication_map.qualifier
    AND annotations.assigned_by = load_go_term_publication_map.assigned_by
    AND annotations.ontology_term_id = load_go_term_publication_map.ontology_term_id
    AND annotations.evidence_code = load_go_term_publication_map.evidence_code
)
ON CONFLICT (go_term_annotation_id, ref_pubmed_id)
DO NOTHING
;
$$,
$$
drop table load_go_term_publication_map;
$$
;
