BEGIN;

-- We want to remove the assemblies we have no regions for. This is because it
-- doesn't make sense to list them in genome browser or anything if we have no
-- data. However, we can't figure out what assemblies to remove until after the
-- release because of the mapping procedure. Even if we don't fetch that data
-- from Ensembl yet we may map sequences into that genome. Thus we wait until
-- after running a release to detect and remove entries from the
-- ensembl_assembly table.
CREATE TEMP TABLE assemblies_with_no_regions as
SELECT
  assem.*
FROM ensembl_assembly AS assem
LEFT JOIN rnc_sequence_regions regions
ON
  regions.assembly_id = assem.assembly_id
WHERE
	regions.id is null
;

DELETE FROM ensembl_assembly
USING assemblies_with_no_regions none
WHERE
  none.assembly_id = ensembl_assembly.assembly_id
;

-- We also want to ensure that all assemblies have example locations. These are
-- used in the genome browser page so they are helpful to have around. We
-- determine interesting ones by hand for important organisms, but otherwise we
-- can just pick one sequence at random. This tries to pick a single sequence.
UPDATE ensembl_assembly assem
SET
  example_chromosome = t.example_chromosome,
  example_start = t.example_start,
  example_end = t.example_end
FROM (
SELECT DISTINCT ON (assem.assembly_id)
	assem.assembly_id,
	regions.chromosome as example_chromosome,
	regions.region_start as example_start,
	regions.region_stop as example_end
FROM ensembl_assembly assem
JOIN rnc_sequence_regions regions
ON
	regions.assembly_id = assem.assembly_id
JOIN rnc_rna_precomputed pre
ON
	pre.id = regions.urs_taxid
WHERE
	assem.example_chromosome is NULL
	AND pre.is_active = true
	AND pre.taxid IS NOT NULL -- Help the query planner a little
ORDER BY assem.assembly_id, regions.chromosome DESC, abs((regions.region_stop - regions.region_start) - 5000) ASC
) t
WHERE
	t.assembly_id = assem.assembly_id
	AND assem.example_chromosome IS NULL
;

COMMIT;
