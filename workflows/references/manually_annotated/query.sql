select
	xref.upi || '_' || xref.taxid,
	-- acc.accession,
	acc."database",
	refs.pmid,
	refs.doi,
	refs.pmcid
	-- refs.epmcid
from rnc_accessions acc
join xref
on xref.ac = acc.accession
join rnc_reference_map rmap on rmap.accession = acc.accession
join rnc_references refs on refs.id = rmap.reference_id
where
	xref.dbid in (24, 20, 14, 16, 18, 23, 27, 44, 48)
	and xref.deleted = 'N'
;
