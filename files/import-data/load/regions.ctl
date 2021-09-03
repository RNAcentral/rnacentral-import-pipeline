LOAD CSV
FROM ALL FILENAMES MATCHING ~<regions.*csv$>
HAVING FIELDS (
  accession,
  region_name,
  chromosome,
  strand,
  assembly_id,
  exon_count,
  exon_start,
  exon_stop
)
INTO {{PGDATABASE}}?load_rnc_sequence_regions
TARGET COLUMNS (
  accession,
  region_name,
  chromosome,
  strand,
  assembly_id,
  exon_count,
  exon_start,
  exon_stop
)

WITH truncate,
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','
;
