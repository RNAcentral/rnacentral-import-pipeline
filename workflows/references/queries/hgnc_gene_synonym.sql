-- HGNC
select
    gene_synonym,
	  gene,
    upi,
    taxid
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" = 'HGNC'
;
