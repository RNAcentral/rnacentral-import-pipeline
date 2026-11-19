-- REFSEQ_CODING: transcript ids and locus tags only; gene symbols would pull in the
-- protein literature. A UTR's id is its mRNA's plus a suffix, so only mRNAs.
select
    urs,
    taxid,
    external_id,
    locus_tag
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" = 'REFSEQ_CODING'
	and ra.feature_name = 'mRNA'
;
