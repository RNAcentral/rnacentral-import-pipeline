COPY (
select
  json_build_object(
    'external_id', acc.external_id,
    'rna_id', xref.upi || '_' || xref.taxid
  )
from xref
join rnc_accessions acc on acc.accession = xref.ac
where
  xref.deleted = 'N'
  and acc.database = 'POMBASE'
) TO STDOUT
