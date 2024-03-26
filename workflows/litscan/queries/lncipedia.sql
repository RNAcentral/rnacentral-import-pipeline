-- LNCipedia (Extract BCRP3 gene)
select
    upi,
    taxid,
    gene,
    gene_synonym -- Split on ,
from xref x
join rnc_accessions ra 
on 
	ra.accession = x.ac 
where
	x.deleted = 'N'
	and ra."database" = 'LNCIPEDIA'
	and (ra."gene" = 'lnc-GGT1-7' or ra."gene" = 'lnc-GGT1-16')
	and ra."gene_synonym" != ''
;
