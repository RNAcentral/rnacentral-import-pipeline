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
        'database', db.descr
    )
FROM xref
join rnc_accessions acc on acc.accession = xref.ac
join rnc_database db on db.id = xref.dbid
join rnc_rna_precomputed pre
on
    pre.upi = xref.upi
    and pre.taxid = xref.taxid
where
    xref.deleted = 'N'
) TO STDOUT
