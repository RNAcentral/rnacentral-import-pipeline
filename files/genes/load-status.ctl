LOAD CSV
FROM ALL FILENAMES MATCHING ~<{{STATUS}}.*csv$>
HAVING FIELDS (
  assembly_id,
  region_id,
  urs_taxid,
  status
)
INTO {{PGDATABASE}}?load_gene_status
TARGET COLUMNS (
  assembly_id,
  region_id,
  urs_taxid,
  status
)

WITH batch rows = 500,
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
