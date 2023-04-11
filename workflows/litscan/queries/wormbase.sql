-- WORMBASE
select
    upi,
    taxid,
    external_id,
    optional_id,
    locus_tag
from xref x
join rnc_accessions ra
on
    ra.accession = x.ac
where
    x.deleted = 'N'
    and ra."database" in ('WORMBASE')
;
