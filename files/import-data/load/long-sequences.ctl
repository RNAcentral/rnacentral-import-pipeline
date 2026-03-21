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
    batch rows = 10000,
    batch size =  64MB,
    prefetch rows = 10000,
    workers = 2,
    concurrency = 1,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '512 MB',
    maintenance_work_mem to '1 GB'
;
