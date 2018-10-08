LOAD CSV
FROM ALL FILENAMES MATCHING ~<secondary_structure.*csv$>
HAVING FIELDS (
    rnc_accession_id,
    secondary_structure,
    md5
) INTO {{PGDATABASE}}
TARGET COLUMNS (
    rnc_accession_id,
    secondary_structure,
    md5
)

WITH
    drop indexes,
    SKIP HEADER = 1,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
drop table if exists load_rnc_secondary_structure;
$$,
$$
create table load_rnc_secondary_structure (
    rnc_accession_id varchar(100),
    secondary_structure text,
    md5 varchar(32)
);
$$
