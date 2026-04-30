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

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_assemblies SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_assemblies;
$$
;
