-- ENA_CODING: the external id is the record accession for a virus and the locus tag for
-- a genome CDS; gene symbols would pull in the protein literature. Only mRNAs, as a
-- UTR's id is its mRNA's plus a suffix.
select
    urs,
    taxid,
    external_id
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" = 'ENA_CODING'
	and ra.feature_name = 'mRNA'
;
