-- GENECARDS
select
    upi,
    taxid,
    gene
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" in ('GENECARDS')
;
