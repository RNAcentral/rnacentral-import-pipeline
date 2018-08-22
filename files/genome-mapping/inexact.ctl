LOAD CSV
FROM '{{FILENAME}}'
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
WHERE identity < 1
AND assembly_id = '{assembly_id}';
$$,
$$ insert into rnc_genome_mapping (rna_id, upi, taxid, region_id, chromosome, "start", stop, strand, assembly_id, "identity") (
    select t1.rna_id, t1.upi, t1.taxid, t1.region_id, t1.chromosome, t1."start", t1.stop, t1.strand, t1.assembly_id, t1."identity"
    from load_genome_mapping t1
    left join rnc_genome_mapping t2
    on t1.rna_id = t2.rna_id
    where t2.rna_id is null
);
; $$,
$$ TRUNCATE TABLE {table}; $$,
$$
UPDATE rnc_rna_precomputed t1
SET has_coordinates = True
FROM rnc_genome_mapping t2
WHERE t2.assembly_id = '{assembly_id}'
AND t1.upi = t2.upi
AND t1.taxid = t2.taxid
; $$
;
"""

