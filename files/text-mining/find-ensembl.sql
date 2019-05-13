COPY (
select
  acc.accession
from rnc_accessions acc
join xref on xref.ac = acc.accession
where
  xref.deleted = 'N'
  and xref.dbid = 25
UNION
select
  acc.external_id
from rnc_accessions acc
join xref on xref.ac = acc.accession
where
  xref.deleted = 'N'
  and xref.dbid = 25
UNION
select
  acc.optional_id
from rnc_accessions acc
join xref on xref.ac = acc.accession
where
  xref.deleted = 'N'
  and xref.dbid = 25
UNION
select
  split_part(acc.optional_id, '.', 1)
from rnc_accessions acc
join xref on xref.ac = acc.accession
where
  xref.deleted = 'N'
  and xref.dbid = 25
) TO STDOUT
