-- ENSEMBL_GENCODE
select
    gene, -- Also search for everything up to the first '.'
    external_id,
    upi,
    taxid
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" in ('ENSEMBL_GENCODE')
;
