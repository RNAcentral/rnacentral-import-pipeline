COPY (
SELECT
  ensembl_url,
  assembly_id,
  taxid,
  division
FROM :species_to_map to_map
) TO STDOUT CSV;
