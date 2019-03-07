COPY (
select
    optional_id
from rnc_accessions acc
join xref on xref.ac = acc.accession
where
    xref.deleted = 'N'
    and xref.dbid = 4
    and optional_id != ''
    and optional_id is not null
) TO STDOUT
