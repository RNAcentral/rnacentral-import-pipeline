COPY (
SELECT
  species.ensembl_url,
  species.assembly_id
FROM ensembl_assembly species
LEFT JOIN rnc_taxonomy tax ON tax.id = species.taxid
WHERE
  exists(select 1 from rnc_sequence_regions reg where reg.assembly_id = species.assembly_id)
  and not coalesce(tax.lineage like '%; Fungi; %', false)
  and species.selected_genome = true
) TO STDOUT CSV;
