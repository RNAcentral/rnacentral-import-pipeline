COPY (
select 
  refs.*
from xref 
join rnc_accessions acc on acc.accession = xref.ac
join rnc_reference_map rmap on rmap.accession = acc.accession
join rnc_references refs on refs.id = rmap.reference_id
where
  xref.deleted = 'N'
  and refs.pmid is not null
) TO STDOUT
