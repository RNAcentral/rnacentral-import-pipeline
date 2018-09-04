LOAD CSV
FROM ALL FILENAMES MATCHING ~<complete_features*.csv>
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
INSERT INTO sequence_features (
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata
) (
SELECT
    upi,
    taxid,
    accession,
    start,
    stop,
    feature_name,
    metadata
from load_crs_features
)
ON CONFLICT (upi, taxid, accession, start, stop, feature_name) DO UPDATE
SET
    metadata = EXCLUDED.metadata
;
$$,
$$
drop table load_crs_features;
$$
;
