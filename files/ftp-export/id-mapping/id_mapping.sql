SET work_mem = '512MB';
SET hash_mem_multiplier = 2;
SET enable_nestloop = off;
SET max_parallel_workers_per_gather = 0;

COPY (
SELECT
    json_build_object(
        'upi', xref.upi,
        'accession', xref.ac,
        'taxid', xref.taxid,
        'external_id', acc.external_id,
        'optional_id', acc.optional_id,
        'rna_type', pre.rna_type,
        'gene', acc.gene,
        'product', acc.product,
        'database', db.descr
    )
FROM xref
join rnc_accessions acc on acc.accession = xref.ac
join rnc_database db on db.id = xref.dbid
join rnc_rna_precomputed pre on pre.id = xref.urs_taxid
where
    xref.deleted = 'N'
    AND right(xref.upi, 1) = ANY (string_to_array(:'chunk', ','))
) TO STDOUT
