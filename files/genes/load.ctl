LOAD CSV
FROM ALL FILENAMES MATCHING ~<locus.*csv$>
HAVING FIELDS (
  assembly_id,
  locus_name,
  locus_public_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  member_count,
  urs_taxid,
  region_id,
  member_status
)
INTO {{PGDATABASE}}?load_locus
TARGET COLUMNS (
  assembly_id,
  locus_name,
  locus_public_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  member_count,
  urs_taxid,
  region_id,
  member_status
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
