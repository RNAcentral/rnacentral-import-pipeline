-- ENSEMBL METAZOA
select
    upi,
    taxid,
    external_id,
    gene, -- Also search for everything up to the first '.'
    locus_tag
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" in ('ENSEMBL_METAZOA')
;
