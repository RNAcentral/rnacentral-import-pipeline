LOAD CSV
FROM ALL FILENAMES MATCHING ~<sequence_regions.*csv$>
HAVING FIELDS (
  accession,
  region_name,
  chromosome,
  strand,
  assembly_id,
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
  exon_start,
  exon_stop
)

WITH truncate,
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rnc_sequence_regions;
$$,
$$
create table load_rnc_sequence_regions (
    accession text,
    urs_taxid text,
    region_name text not null,
    chromosome text,
    strand int4,
    exon_start int4,
    exon_stop int4,
    assembly_id varchar(255),
    providing_database text
);
$$
;
