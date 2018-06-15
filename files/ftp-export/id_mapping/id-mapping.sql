SELECT
    xref.upi,
    xref.ac AS accession,
    xref.taxid,
    acc.external_id,
    acc.optional_id,
    pre.rna_type,
    acc.gene,
    db.descr AS database
FROM xref
join rnc_accessions acc on acc.accession = xref.ac
join rnc_database db on db.id = xref.dbid
join rnc_rna_precomputed pre on pre.upi = xref.upi and pre.taxid = xref.taxid
where
    xref.deleted = 'N'

