-- SGD_CODING: systematic ORF names, e.g. YAL001C. UTRs share them, so only mRNAs.
select
    urs,
    taxid,
    locus_tag
from xref x
join rnc_accessions ra
on
	ra.accession = x.ac
where
	x.deleted = 'N'
	and ra."database" = 'SGD_CODING'
	and ra.feature_name = 'mRNA'
;
