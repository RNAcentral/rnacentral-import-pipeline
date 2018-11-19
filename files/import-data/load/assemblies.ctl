LOAD CSV
FROM ALL FILENAMES MATCHING ~<assemblies.*csv>
HAVING FIELDS (
    assembly_id,
    assembly_full_name,
    gca_accession,
    assembly_ucsc,
    common_name,
    taxid,
    ensembl_url,
    division,
    subdomain,
    example_chromosome,
    example_start,
    example_end,
    blat_mapping
)
INTO {{PGDATABASE}}?load_assemblies
TARGET COLUMNS (
    assembly_id,
    assembly_full_name,
    gca_accession,
    assembly_ucsc,
    common_name,
    taxid,
    ensembl_url,
    division,
    subdomain,
    example_chromosome,
    example_start,
    example_end,
    blat_mapping
)

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_assemblies;
$$,
$$
CREATE TABLE load_assemblies (
	assembly_id varchar(255) NOT NULL,
	assembly_full_name varchar(255) NOT NULL,
	gca_accession varchar(20) NULL,
	assembly_ucsc varchar(100) NULL,
	common_name varchar(255),
	taxid int4 NOT NULL,
	ensembl_url varchar(100) NULL,
	division varchar(20) NULL,
	blat_mapping int4 NULL,
	example_chromosome varchar(20) NULL,
	example_end int4 NULL,
	example_start int4 NULL,
	subdomain varchar(100) NOT NULL
);
$$
;
