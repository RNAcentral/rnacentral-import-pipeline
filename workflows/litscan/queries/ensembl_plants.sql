-- ENSEMBL PLANTS
select
    upi,
    taxid,
    external_id,
    gene,
    locus_tag
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" in ('ENSEMBL_PLANTS')
;
