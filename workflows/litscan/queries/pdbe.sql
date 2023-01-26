-- PDBE
select
    upi,
    taxid,
	  external_id
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" = 'PDBE'
;
