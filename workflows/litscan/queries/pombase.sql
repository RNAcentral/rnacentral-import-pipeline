-- POMBASE
select
    upi,
    taxid,
    external_id,
	  gene,
	  gene_synonym -- SPlit on ','
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" in ('POMBASE')
;
