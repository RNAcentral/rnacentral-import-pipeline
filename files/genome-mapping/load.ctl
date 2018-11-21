LOAD CSV
FROM ALL FILENAMES MATCHING ~<genome-mapping.*csv$>
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
INTO {{DB_URL}}?load_genome_mapping
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

WITH
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_genome_mapping;
$$,
$$
create table load_genome_mapping (
    accession text,
    urs_taxid text,
    region_name text not null,
    chromosome text,
    strand int4,
    exon_start int4,
    exon_stop int4,
    assembly_id varchar(255),
    exon_count int,
    identity float,
    providing_database text
);
$$

AFTER LOAD DO
$$
DELETE FROM rnc_sequence_regions regions
USING load_genome_mapping load
WHERE
    regions.assembly_id = load.assembly_id
    AND regions.was_mapped = true
;
$$,
$$
INSERT INTO rnc_sequence_regions (
    urs_taxid,
    region_name,
    chromosome,
    strand,
    region_start,
    region_stop,
    assembly_id,
    exon_count,
    was_mapped,
    identity,
    providing_databases
) (
SELECT
    max(load.urs_taxid),
    load.region_name,
    max(load.chromosome),
    max(load.strand),
    min(load.exon_start),
    max(load.exon_stop),
    load.assembly_id,
    max(load.exon_count),
    true,
    load.identity,
    '{}'::text[]
FROM load_genome_mapping load
JOIN ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
GROUP BY load.region_name, load.assembly_id
) ON CONFLICT (region_name, assembly_id) DO NOTHING
;
$$,
$$
INSERT INTO rnc_sequence_exons (
    region_id,
    exon_start,
    exon_stop
) (
SELECT
    regions.id,
    load.exon_start,
    load.exon_stop
FROM load_genome_mapping load
JOIN rnc_sequence_regions regions ON regions.region_name = load.region_name
JOIN ensembl_assembly ensembl ON ensembl.assembly_id = load.assembly_id
) ON CONFLICT (region_id, exon_start, exon_stop) DO NOTHING
;
$$,
$$
DROP TABLE load_genome_mapping;
$$
;
