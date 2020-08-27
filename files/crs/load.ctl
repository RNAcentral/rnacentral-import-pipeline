LOAD CSV
FROM ALL FILENAMES MATCHING ~<complete_features.*csv$>
HAVING FIELDS (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata
)
INTO {{PGDATABASE}}?load_crs_features
TARGET COLUMNS (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata
)
WITH
    drop indexes,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_crs_features;
$$,
$$
CREATE TABLE load_crs_features (
    upi text,
    taxid int,
    accession text,
    start int,
    stop int,
    feature_name text,
    metadata jsonb
);
$$

AFTER LOAD DO
$$
delete from rnc_sequence_features where feature_name = 'conserved_rna_structure';
$$,
$$
INSERT INTO rnc_sequence_features (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata
) (
SELECT distinct
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata
from load_crs_features
)
;
$$,
$$
drop table load_crs_features;
$$
;

