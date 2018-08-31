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
INTO {{PGDATABASE}}?rnc_sequence_features
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

SET
    work_mem to '256 MB',
    maintenance_work_mem to '256 GB',
;
