COPY (
SELECT
	DISTINCT regions.assembly_id, lower(db.descr)
FROM rnc_database db
JOIN xref ON xref.dbid = db.id
JOIN rnc_sequence_regions regions
ON
  regions.urs_taxid = xref.upi || '_' || xref.taxid
JOIN ensembl_coordinate_systems coords
ON
  coords.chromosome = regions.chromosome
  AND coords.assembly_id = regions.assembly_id
WHERE
	db.descr IN ('FLYBASE', 'POMBASE', 'TAIR')
  AND coords.is_reference = true
) TO STDOUT CSV

-- This uses the regions and not a join on xref.taxid to the assembly table
-- because I think it will be more reliable to use the regions. Sometimes we
-- have to ust a different 'correct' the taxid to properly connect sequences to
-- organisms
