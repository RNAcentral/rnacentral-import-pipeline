LOAD CSV
FROM ALL FILENAMES MATCHING ~<long_sequences.*csv$>
HAVING FIELDS (
    CRC64,
    LEN,
    seq_long,
    DATABASE,
    AC,
    OPTIONAL_ID,
    VERSION,
    TAXID,
    MD5
)
INTO {{PGDATABASE}}?load_rnacentral_all
TARGET COLUMNS (
    CRC64,
    LEN,
    seq_long,
    DATABASE,
    AC,
    OPTIONAL_ID,
    VERSION,
    TAXID,
    MD5
)

WITH
    drop indexes,
    batch rows = 25000,
    batch size =  512MB,
    workers = 10,
    concurrency = 2,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '256 MB',
    maintenance_work_mem to '1 GB'
;
