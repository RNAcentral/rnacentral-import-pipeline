COPY (
select
	distinct regions.assembly_id, lower(db.descr)
from rnc_database db
join xref on xref.dbid = db.id
join rnc_sequence_regions regions on regions.urs_taxid = xref.upi || '_' || xref.taxid
where
	db.descr in ('FLYBASE', 'POMBASE', 'TAIR')
) TO STDOUT CSV

-- This uses the regions and not a join on xref.taxid to the assembly table
-- because I think it will be more reliable to use the regions. Sometimes we
-- have to ust a different 'correct' the taxid to properly connect sequences to
-- organisms
