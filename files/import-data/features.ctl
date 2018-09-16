LOAD CSV
FROM ALL FILENAMES MATCHING ~<features.*csv$>
HAVING FIELDS (
    accession,
    taxid,
    start,
    stop,
    feature_name,
    metadata
)
INTO {{PGDATABASE}}?load_rnc_sequence_features
TARGET COLUMNS (
    accession,
    taxid,
    start,
    stop,
    feature_name,
    metadata
)

WITH truncate,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rnc_sequence_features;
$$,
$$
create table load_rnc_sequence_features (
    accession varchar(100) NOT NULL,
    taxid int not null,
    start int not null,
    stop int not null,
    feature_name varchar(50),
    metadata jsonb
);
$$
;
