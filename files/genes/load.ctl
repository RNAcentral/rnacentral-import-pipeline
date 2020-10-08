LOAD CSV
FROM ALL FILENAMES MATCHING ~<genes.*csv$>
HAVING FIELDS (
  taxid,
  assembly_id,
  locus_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  urs_taxid,
  region_id,
  is_representative
)
INTO {{PGDATABASE}}?load_locus
TARGET COLUMNS (
  taxid,
  assembly_id,
  locus_name,
  chromosome,
  strand text,
  locus_start,
  locus_stop,
  urs_taxid,
  region_id,
  is_representative
)

WITH truncate,
    drop indexes,
    batch rows = 500,
    batch size = 32MB,
    prefetch rows = 500,
    workers = 2, concurrency = 1,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '256 MB',
    maintenance_work_mem to '256 GB'
;
