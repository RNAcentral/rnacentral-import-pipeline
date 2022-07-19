-- REFSEQ
select
    upi,
    taxid,
    external_id,
    gene,
    gene_synonym, -- Split on ','
    optional_id
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" in ('REFSEQ')
;
