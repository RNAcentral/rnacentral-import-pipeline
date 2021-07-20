LOAD CSV
FROM ALL FILENAMES MATCHING ~<ensembl-pseudogenes.*csv$>
HAVING FIELDS (
  gene,
  region_name,
  chromosome,
  strand,
  assembly_id,
  exon_count,
  exon_start,
  exon_stop
)
INTO {{PGDATABASE}}?load_ensembl_pseudogenes
TARGET COLUMNS (
  gene,
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
