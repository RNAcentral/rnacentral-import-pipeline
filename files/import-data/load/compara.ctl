LOAD CSV
FROM ALL FILENAMES MATCHING ~<compara.*csv>
HAVING FIELDS (
    homology_group,
    ensembl_transcript
)
INTO {{PGDATABASE}}?load_compara
TARGET COLUMNS (
    homology_group,
    ensembl_transcript
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_compara;
$$,
$$
CREATE TABLE load_compara (
	homology_group text not null,
  ensembl_transcript text not null
);
$$
;
