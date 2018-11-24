LOAD CSV
FROM ALL FILENAMES MATCHING ~<karyotypes.*csv>
HAVING FIELDS (
    assembly_id,
    karyotype
)
INTO {{PGDATABASE}}?load_karyotypes
TARGET COLUMNS (
    assembly_id,
    karyotype
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_karyotypes;
$$,
$$
CREATE TABLE load_karyotypes (
	assembly_id varchar(255) NOT NULL,
    karyotype text
);
$$
;
