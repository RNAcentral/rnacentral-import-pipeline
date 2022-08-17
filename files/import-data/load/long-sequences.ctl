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
    batch size =  256MB,
    prefetch rows = 50000,
    workers = 5,
    concurrency = 2,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '512 MB',
    maintenance_work_mem to '1 GB'
;
