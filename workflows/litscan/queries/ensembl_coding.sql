-- ENSEMBL_CODING: transcript and gene ids only; protein-coding gene symbols would
-- pull in the protein literature. UTRs share their mRNA's ids, so only mRNAs.
select
    urs,
    taxid,
    external_id,
    gene
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" = 'ENSEMBL_CODING'
	and ra.feature_name = 'mRNA'
;
