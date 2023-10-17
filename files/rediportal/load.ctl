LOAD CSV
FROM ALL FILENAMES MATCHING ~<rediportal-data.*csv$>
HAVING FIELDS (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata,
    feature_provider
)
INTO {{PGDATABASE}}?load_rediportal_features
TARGET COLUMNS (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata,
    feature_provider
)
WITH
    drop indexes,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rediportal_features;
$$,
$$
CREATE TABLE load_rediportal_features (
    upi text,
    taxid int,
    accession text,
    start int,
    stop int,
    feature_name text,
    metadata jsonb,
    feature_provider text
);
$$

AFTER LOAD DO
$$
delete from rnc_sequence_features where feature_name = 'rna_editing_event';
$$,
$$
INSERT INTO rnc_sequence_features (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata,
    feature_provider
) (
SELECT distinct
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata,
    feature_provider
from load_rediportal_features
)
;
$$,
$$
drop table load_rediportal_features;
$$
;
