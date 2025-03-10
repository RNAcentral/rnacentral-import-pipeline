LOAD CSV
FROM ALL FILENAMES MATCHING ~<accessions.*csv$>
HAVING FIELDS (
    accession,
    parent_ac,
    seq_version,
    feature_start,
    feature_end,
    feature_name,
    non_coding_id,
    database,
    external_id,
    optional_id,
    description,
    organelle,
    chromosome,
    function,
    gene,
    gene_synonym,
    inference,
    locus_tag,
    mol_type,
    ncrna_class,
    note,
    product,
    standard_name,
    db_xref,
    so_term,
    url
)
INTO {{PGDATABASE}}?load_rnc_accessions
TARGET COLUMNS (
    accession,
    parent_ac,
    seq_version,
    feature_start,
    feature_end,
    feature_name,
    non_coding_id,
    database,
    external_id,
    optional_id,
    description,
    organelle,
    chromosome,
    function,
    gene,
    gene_synonym,
    inference,
    locus_tag,
    mol_type,
    ncrna_class,
    note,
    product,
    standard_name,
    db_xref,
    so_term,
    url
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

BEFORE LOAD DO
$$
truncate table load_rnc_accessions;
$$,
$$
ALTER TABLE rnacen.load_rnc_accessions SET (
    autovacuum_enabled = false,
    toast.autovacuum_enabled = false
);
$$

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_accessions SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$
;
