LOAD CSV
FROM '{{PATTERN}}'
HAVING FIELDS
(
    rna_id,
    upi,
    taxid,
    unused_field,
    region_id,
    chromosome,
    start,
    stop,
    strand,
    assembly_id,
    identity
)
INTO {{DB_URL}}
TARGET COLUMNS
(
    rna_id,
    upi,
    taxid,
    region_id,
    chromosome,
    start,
    stop,
    strand,
    assembly_id,
    identity
)
WITH
    fields terminated by '\t'
SET
    search_path = 'rnacen'
BEFORE LOAD DO
$$ DROP TABLE IF EXISTS {table}; $$,
$$ CREATE TABLE IF NOT EXISTS rnacen.{table} (
    rna_id varchar(50) NULL,
    upi varchar(13) NULL,
    taxid int4 NULL,
    region_id varchar(100) NULL,
    chromosome varchar(100) NULL,
    "start" int4 NULL,
    stop int4 NULL,
    strand varchar(10) NULL,
    assembly_id varchar(50) NULL,
    identity float8
); $$
AFTER LOAD DO
$$ DELETE FROM rnc_genome_mapping
WHERE identity = 1
AND assembly_id IN (SELECT distinct assembly_id from {table});
$$,
$$ INSERT INTO rnc_genome_mapping (rna_id, upi, taxid, region_id, chromosome,
                                   "start", stop, strand, assembly_id, "identity")
    (SELECT rna_id, upi, taxid, region_id, chromosome, "start", stop, strand,
            assembly_id, "identity"
	from {table})
; $$
;
